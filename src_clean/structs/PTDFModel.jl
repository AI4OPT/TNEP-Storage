using JuMP
using CSV, DataFrames
include("OptimizationModel.jl")
include("../helpers/helpers.jl")

# PTDF Model with iterative constraint addition
mutable struct PTDFModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String
    
    # PTDF-specific fields
    ptdf_matrix::Matrix{Float64}
    tracked_constraints::Dict{Tuple{Int,Int,Int,Bool}, Bool}
    rate_a_nonzero::Set{Int}
    max_ptdf_iterations::Int
    max_ptdf_per_iteration::Int
    ptdf_tol::Float64
    solve_time::Float64
    
    # Pre-computed mappings for efficiency
    gen_bus_map::Dict{Int, Int}
    
    # Dimension sizes for convenience
    R::Int  # num_representatives
    N::Int  # num buses
    E::Int  # num branches
    T::Int  # num hours
    G::Int  # num generators
    
    # Inner constructor
    function PTDFModel(data::Dict{String, Any}, optimizer, simdir::String; 
                      max_ptdf_iterations::Int=256,
                      max_ptdf_per_iteration::Int=32,
                      ptdf_tol::Float64=1e-6)
        
        # Create the base JuMP model
        jump_model = create_base_model(data, optimizer)
        
        # Compute PTDF matrix
        ptdf_matrix = do_all_ptdf(data)
        
        # Apply cutoff if specified
        if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
            ptdf_sparse = map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"]
            ptdf_matrix = ptdf_matrix .* ptdf_sparse
        end
        
        # Get rate_a nonzero branches
        rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
        rate_a_set = Set(parse(Int, x) for x in rate_a_nonzero)
        
        # Pre-compute gen_bus_map
        gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] 
                          for g in keys(data["gen"]))
        
        # Get dimension sizes
        R = data["param"]["num_representatives"]
        N = length(data["bus"])
        E = length(data["branch"])
        T = data["param"]["num_hours"]
        G = length(data["gen"])
        
        # Add container for flow constraints to JuMP model
        jump_model[:ptdf_flow] = Dict{String, ConstraintRef}()
        
        # Create and return the model
        new(jump_model, data, simdir, ptdf_matrix, 
            Dict{Tuple{Int,Int,Int,Bool}, Bool}(),
            rate_a_set, max_ptdf_iterations, max_ptdf_per_iteration, 
            ptdf_tol, 0.0, gen_bus_map, R, N, E, T, G)
    end
end

# Creates skeleton PTDF model
function create_base_model(data::Dict{String, Any}, optimizer)    
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    model = JuMP.Model(optimizer)

    #
    #   I. Variables
    #
    
    # Generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    @variable(model, pg[r in 1:R, g in 1:G, t in 1:T])

    for g in 1:G
        gen = data["gen"]["$g"]
        is_renewable = gen["gen_type"] ∈ data["param"]["renewable_types"]
        is_foreign = gen["gen_type"] ∈ ["foreign"]
        if is_renewable
            set_lower_bound.(pg[:, g, :], 0.0)
            for r in 1:R, t in 1:T
                set_upper_bound(pg[r, g, t], max(0, gen["profile"]["$r"][t]))
            end
        elseif is_foreign
            for r in 1:R, t in 1:T
                fix(pg[r, g, t], gen["profile"]["$r"][t])
            end
        else
            set_lower_bound.(pg[:, g, :], gen["pmin"])
            set_upper_bound.(pg[:, g, :], gen["pmax"])
        end
    end

    # Conditional variable declaration based on whether relaxed model is used
    if haskey(data["param"], "relaxed_first_stage") && data["param"]["relaxed_first_stage"] == true
        # Use continuous variables for relaxed model
        @variable(model, 0 <= gamma[a=1:E] <= K)
    else
        # Use integer variables for standard model
        @variable(model, 0 <= gamma[a=1:E] <= K, Int)
    end
    
    @variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0)  # under-served energy at bus
    @variable(model, s_energy[i=1:N] >= 0)  # energy rating of storage
    @variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0)  # state of charge of storage
    @variable(model, ch[r=1:R, i=1:N, t=1:T] >= 0)  # charging of storage
    @variable(model, dis[r=1:R, i=1:N, t=1:T] >= 0)  # discharging of storage

    #
    #   II. Constraints
    #
    
    # Global power balance
    @constraint(model, 
        power_balance[r in 1:R, t in 1:T],
        sum(pg[r,g,t] for g in 1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i in 1:N)
        + sum(ue[r,i,t] for i in 1:N) 
        - sum(ch[r,i,t] for i in 1:N)
        + sum(dis[r,i,t] for i in 1:N) 
        == 0
    )

    # Check if previous investments exist
    if haskey(data["param"], "previous_investment_dir")
        prev_dir = data["param"]["previous_investment_dir"]
        add_prev_upgrades(model, data, gamma, s_energy, prev_dir)
    end

    # Check if current investments to test
    if haskey(data["param"], "current_investment_dir")
        cur_dir = data["param"]["current_investment_dir"]
        trans_file = joinpath(cur_dir, "line_investments.csv")
        storage_file = joinpath(cur_dir, "storage_investments.csv")
        
        # Only add transmission constraints if the file exists
        if isfile(trans_file)
            trans_df = CSV.read(trans_file, DataFrame)
            @constraint(model,
                current_gamma[a in 1:E],
                gamma[a] == trans_df[a, :Upgrade_Lvl]
            )
        end
        
        # Only add storage constraints if the file exists
        if isfile(storage_file)
            storage_df = CSV.read(storage_file, DataFrame)
            if get(data["param"], "storage_needs_scaling", false)
                @constraint(model,
                    current_s_energy[i in 1:N],
                    s_energy[i] == storage_df[i, :Storage_Energy] * data["param"]["storage_energy_size"]
                )
            else
                @constraint(model,
                    current_s_energy[i in 1:N],
                    s_energy[i] == storage_df[i, :Storage_Energy]
                )
            end
        end
    end

    # Check if fixed investments to test (alternative parameter name)
    if haskey(data["param"], "inv_dir")
        inv_dir = data["param"]["inv_dir"]
        add_prev_upgrades(model, data, gamma, s_energy, inv_dir; equality=true)
    end

    # If rate a is zero (unlimited), then don't allow upgrades
    @constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # SOC over time constraint
    @constraint(model, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - 
                      dis[r,i,t] / data["param"]["bess_efficiency"]
    )

    # SOC init and end constraint
    @constraint(model,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] + 
                      ch[r,i,1] * data["param"]["bess_efficiency"] - 
                      dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    @constraint(model,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i]
    )

    # SOC energy rating constraint
    @constraint(model, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i]
    )

    # Energy rating only if storage installed
    @constraint(model, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= data["param"]["max_energy_rating"]
    )

    # Charge/discharge must be constrained by power rating
    @constraint(model,
        charge_discharge_lb[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] <= s_energy[i] / 4
    )
    @constraint(model,
        charge_discharge_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= s_energy[i] / 4
    )

    # OPTIONAL: precompute maximum nodal injections
    if get(data["param"], "nodecap", false)
        nodal_caps = Dict{Int, Float64}()
        for i in 1:N
            max_outflow = sum(data["branch"]["$a"]["rate_a"] + 
                            get_capacity_increment(data, a) * K 
                            for a in data["arcs_from"]["$i"])
            
            max_charge = min(data["param"]["max_power_rating"], 
                           data["param"]["max_energy_rating"] / 4)
            nodal_caps[i] = max_outflow + max_charge                
        end

        @constraint(model,
            nodal_injection_cap[r in 1:R, i in 1:N, t in 1:T],
            sum(model[:pg][r, g, t] for g in 1:G if gen_bus_map[g] == i) <=
            data["bus"]["$i"]["load"]["$r"][t] + nodal_caps[i]
        )
    end
    
    #
    #   III. Objective
    #
    
    operational_weight = get(data["param"], "operational_weight", 1)

    objective_expr = (
        sum(s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
        sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * 
            data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i in 1:N)
                for t in 1:T)
            )
        for r in 1:R) * operational_weight
    )

    if haskey(data["param"], "only_feasibility") && data["param"]["only_feasibility"]
        @objective(model, Min, objective_expr)
    else
        objective_expr += sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g in 1:G) +
                    sum(get(data["param"], "storage_operation_cost", 0.0) * 
                        (ch[r,i,t] + dis[r,i,t]) for i=1:N)
                for t in 1:T)
            )
        for r in 1:R) * operational_weight
        @objective(model, Min, objective_expr)
    end
    
    return model
end

# Warm-start from saved constraints
function warm_start!(model::PTDFModel)
    constraint_file = joinpath(model.simdir, "tracked_constraints.csv")
    if !isfile(constraint_file)
        return
    end
    
    tracked_df = CSV.read(constraint_file, DataFrame)
    println("Beginning optimized warm-start process")
    
    # Group the constraints by arc once
    arc_groups = groupby(tracked_df, :arc)
    total_constraints = 0
    
    for group in arc_groups
        a = first(group).arc
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        cap_increment = get_capacity_increment(model.data, a)
        ptdf_row = model.ptdf_matrix[a, :]
        
        # Extract all the (r,t) pairs for this arc
        rep_time_pairs = [(cons.rep, cons.time, cons.ub) for cons in eachrow(group)]
        
        # Pre-compute the line limits for all constraints in this group
        line_limit = line_limit_base + model.jump_model[:gamma][a] * cap_increment
        
        # Create constraints for each (r,t) pair
        for (r, t, ub) in rep_time_pairs
            # Build the flow expression
            flow_expr = sum(ptdf_row[i] * (
                sum(model.jump_model[:pg][r, g, t] for g in 1:model.G 
                    if model.gen_bus_map[g] == i; init=0.0) -
                model.data["bus"]["$i"]["load"]["$r"][t] +
                model.jump_model[:ue][r, i, t] -
                model.jump_model[:ch][r, i, t] +
                model.jump_model[:dis][r, i, t]
            ) for i in 1:model.N)
            
            # Add constraint
            if ub
                model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = 
                    @constraint(model.jump_model, flow_expr <= line_limit)
            else
                model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = 
                    @constraint(model.jump_model, flow_expr >= -line_limit)
            end
            
            # Track the constraint
            model.tracked_constraints[(a, r, t, ub)] = true
            total_constraints += 1
        end
    end
    
    println("Warm-started with $total_constraints constraints")
end

# Find violations in current solution
function find_violations(model::PTDFModel, flows, gamma_values)
    violations = []
    
    for a in model.rate_a_nonzero
        # Pre-compute for this branch
        cap_increment = get_capacity_increment(model.data, a)
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment
        
        for r in 1:model.R, t in 1:model.T
            flow_val = flows[a, r, t]
            violation_amount = abs(flow_val) - line_limit_fixed
            ub = (flow_val >= 0)

            # Skip if already tracked
            if get(model.tracked_constraints, (a, r, t, ub), false)
                continue
            end
    
            if violation_amount > model.ptdf_tol
                push!(violations, (a, r, t, ub, violation_amount))
            end
        end
    end
    
    return violations
end

# Add PTDF constraints for violations
function add_constraints!(model::PTDFModel, sorted_violations)
    n_added = 0
    max_to_add = model.max_ptdf_per_iteration
    
    for (a, r, t, ub, _) in Iterators.take(sorted_violations, max_to_add)
        ptdf_row = model.ptdf_matrix[a, :]
        cap_increment = get_capacity_increment(model.data, a)
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        line_limit_expr = line_limit_base + model.jump_model[:gamma][a] * cap_increment
    
        # Define the flow expression
        flow_expr = sum(ptdf_row[i] * (
            sum(model.jump_model[:pg][r, g, t] for g in 1:model.G 
                if model.gen_bus_map[g] == i; init=0.0)
            - model.data["bus"]["$i"]["load"]["$r"][t]
            + model.jump_model[:ue][r, i, t] 
            - model.jump_model[:ch][r, i, t] 
            + model.jump_model[:dis][r, i, t]
        ) for i in 1:model.N)
    
        # Add constraint
        if ub
            model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = 
                @constraint(model.jump_model, flow_expr <= line_limit_expr)
        else
            model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = 
                @constraint(model.jump_model, flow_expr >= -line_limit_expr)
        end
    
        model.tracked_constraints[(a, r, t, ub)] = true
        n_added += 1
    end
    
    return n_added
end

# Main solve function with iterative constraint addition
function solve!(model::PTDFModel; configure_optimizer::Bool=true)
    # Configure optimizer if requested
    if configure_optimizer
        set_optimizer_attribute(model.jump_model, "LogFile", 
                               joinpath(model.simdir, "gurobi_logfile.log"))
        set_optimizer_attribute(model.jump_model, "MIPGap", 
                               model.data["param"]["mip_gap"])
    end
    
    # Warm-start from previous constraints if available
    warm_start!(model)
    
    # Begin lazy PTDF loop
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < model.max_ptdf_iterations
        # Solve model
        JuMP.optimize!(model.jump_model)
        
        # Exit if not solved optimally
        st = termination_status(model.jump_model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(model.jump_model[:pg])
        ue_values = value.(model.jump_model[:ue])
        ch_values = value.(model.jump_model[:ch])
        dis_values = value.(model.jump_model[:dis])
        gamma_values = value.(model.jump_model[:gamma])
        
        # Compute all flows at once
        flows = zeros(length(model.rate_a_nonzero), model.R, model.T)
        compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, 
                      model.data, model.ptdf_matrix)

        # Find violations
        violations = find_violations(model, flows, gamma_values)
        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        
        # Add constraints for top violations
        n_added = add_constraints!(model, sorted_violations)

        if n_added > 0
            save_tracked_constraints(model.simdir, model, n_violated)
        end
        
        solved = (n_violated == 0)
        niter += 1
        println("Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    model.solve_time = time() - t0
    save_solve_time(model.simdir, model.solve_time)
    save_power_injections(model.simdir, model.jump_model, model.data)
    
    return model.jump_model
end

# Helper function to compute flows
function compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, data, PTDF)
    n_buses = length(data["bus"])
    net_injections = zeros(n_buses)
    
    for r in 1:size(pg_values, 1), t in 1:size(pg_values, 3)
        fill!(net_injections, 0.0)

        for i in 1:n_buses 
            net_injections[i] = sum(pg_values[r, g, t] 
                                   for g in 1:length(data["gen"]) 
                                   if data["gen"]["$g"]["gen_bus"] == i; 
                                   init=0.0)
            net_injections[i] -= data["bus"]["$i"]["load"]["$r"][t]
            net_injections[i] += ue_values[r, i, t]
            net_injections[i] -= ch_values[r, i, t]
            net_injections[i] += dis_values[r, i, t]
        end
        
        flows[:, r, t] = PTDF * net_injections
    end
end

function save_tracked_constraints(simdir, model::PTDFModel, n_violated)
    """
    Save tracked constraints from PTDFModel to CSV file.
    """
    # Convert dictionary to DataFrame
    df = DataFrame(
        arc = [k[1] for k in keys(model.tracked_constraints)],
        rep = [k[2] for k in keys(model.tracked_constraints)],
        time = [k[3] for k in keys(model.tracked_constraints)],
        ub = [k[4] for k in keys(model.tracked_constraints)],
        tracked = collect(values(model.tracked_constraints)),
        violations_left = [n_violated for i in keys(model.tracked_constraints)]
    )
    
    # Ensure output directory exists
    output_dir = joinpath(simdir, "output")
    mkpath(output_dir)
    
    # Save to CSV
    CSV.write(joinpath(output_dir, "tracked_constraints.csv"), df)
end

# Overloaded version that accepts PTDFModel directly
function save_power_injections(simdir, model::PTDFModel)
    """
    Save power injections from PTDFModel.
    """
    save_power_injections(simdir, model.jump_model, model.data)
end

# Convenience constructor that matches old function signature
function create_model_r1_ptdf_iterative_simplified_sorted(
    simdir, data::Dict{String, Any}, optimizer; 
    max_ptdf_iterations::Int=256,
    max_ptdf_per_iteration::Int=32,
    ptdf_tol::Float64=1e-6)
    
    model = PTDFModel(data, optimizer, simdir;
                     max_ptdf_iterations=max_ptdf_iterations,
                     max_ptdf_per_iteration=max_ptdf_per_iteration,
                     ptdf_tol=ptdf_tol)    
    return model
end