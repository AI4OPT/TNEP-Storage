using JuMP
using CSV, DataFrames
include("OptimizationModel.jl")
include("../helpers/helpers.jl")

mutable struct ConstraintStat
    first_iter::Int          # iteration where this constraint was first violated
    total_appearances::Int   # number of iterations where it was violated
    max_violation::Float64   # max (abs(flow) - line_limit) observed
    tightness_history::Vector{Tuple{Int, Float64}}  # (iter, abs(flow)-line_limit): neg=slack, pos=violated
end

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
    
    # Violation statistics for pruning analysis
    constraint_stats::Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}

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
                      max_ptdf_per_iteration::Int=-1,
                      ptdf_tol::Float64=1e-6)

        max_ptdf_per_iteration = max_ptdf_per_iteration == -1 ?
            get(data["param"], "max_ptdf_per_iteration", 32) : max_ptdf_per_iteration

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
        model = new(jump_model, data, simdir, ptdf_matrix,
            Dict{Tuple{Int,Int,Int,Bool}, Bool}(),
            rate_a_set, max_ptdf_iterations, max_ptdf_per_iteration,
            ptdf_tol, 0.0,
            Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}(),
            gen_bus_map, R, N, E, T, G)

        if haskey(data["param"], "previous_investment_dir")
            prev_dir = data["param"]["previous_investment_dir"]
            gamma_min, s_energy_min = load_investments_from_dir(
                prev_dir, data, E=E, N=N, allow_missing=false
            )
            set_previous_investments!(model, gamma_min, s_energy_min)
        end
    
        if haskey(data["param"], "current_investment_dir")
            cur_dir = data["param"]["current_investment_dir"]
            gamma_val, s_energy_val = load_investments_from_dir(
                cur_dir, data, E=E, N=N, allow_missing=true
            )
            fix_investments!(model, gamma_val, s_energy_val)
        end
        
        return model
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

    jump_model = JuMP.Model(optimizer)

    #
    #   I. Variables
    #
    
    # Generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    @variable(jump_model, pg[r in 1:R, g in 1:G, t in 1:T])

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
        @variable(jump_model, 0 <= gamma[a=1:E] <= K)
    else
        # Use integer variables for standard model
        @variable(jump_model, 0 <= gamma[a=1:E] <= K, Int)
    end
    
    @variable(jump_model, ue[r=1:R, i=1:N, t=1:T] >= 0)  # under-served energy at bus
    @variable(jump_model, s_energy[i=1:N] >= 0)  # energy rating of storage
    @variable(jump_model, soc[r=1:R, i=1:N, t=1:T] >= 0)  # state of charge of storage
    @variable(jump_model, ch[r=1:R, i=1:N, t=1:T] >= 0)  # charging of storage
    @variable(jump_model, dis[r=1:R, i=1:N, t=1:T] >= 0)  # discharging of storage

    #
    #   II. Constraints
    #
    
    # Global power balance
    @constraint(jump_model, 
        power_balance[r in 1:R, t in 1:T],
        sum(pg[r,g,t] for g in 1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i in 1:N)
        + sum(ue[r,i,t] for i in 1:N) 
        - sum(ch[r,i,t] for i in 1:N)
        + sum(dis[r,i,t] for i in 1:N) 
        == 0
    )

    # If rate a is zero (unlimited), then don't allow upgrades
    @constraint(jump_model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # SOC over time constraint
    @constraint(jump_model, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == (1 - get(data["param"], "self_discharge", 0)) * soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - 
                      dis[r,i,t] / data["param"]["bess_efficiency"]
    )

    # SOC init and end constraint
    @constraint(jump_model,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] * data["param"]["storage_energy_size"] + 
                      ch[r,i,1] * data["param"]["bess_efficiency"] - 
                      dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    @constraint(jump_model,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] * data["param"]["storage_energy_size"]
    )

    # SOC energy rating constraint
    @constraint(jump_model, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i] * data["param"]["storage_energy_size"]
    )

    # Energy rating only if storage installed
    @constraint(jump_model, 
        installed_energy_ub[i in 1:N],
        s_energy[i] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
    )

    # Charge/discharge must be constrained by power rating
    @constraint(jump_model,
        charge_discharge_lb[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] <= s_energy[i] * data["param"]["storage_energy_size"] / 4
    )
    @constraint(jump_model,
        charge_discharge_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= s_energy[i] * data["param"]["storage_energy_size"] / 4
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

        @constraint(jump_model,
            nodal_injection_cap[r in 1:R, i in 1:N, t in 1:T],
            sum(jump_model[:pg][r, g, t] for g in 1:G if gen_bus_map[g] == i) <=
            data["bus"]["$i"]["load"]["$r"][t] + nodal_caps[i]
        )
    end
    
    #
    #   III. Objective
    #
    
    operational_weight = get(data["param"], "operational_weight", 1)

    objective_expr = (
        sum(s_energy[i] * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"] for i in 1:N) +
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
        @objective(jump_model, Min, objective_expr)
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
        @objective(jump_model, Min, objective_expr)
    end
    
    return jump_model
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

# Find violations in current solution; niter is 1-indexed current iteration
function find_violations(model::PTDFModel, flows, gamma_values, niter::Int)
    violations = []

    for a in model.rate_a_nonzero
        cap_increment = get_capacity_increment(model.data, a)
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment

        for r in 1:model.R, t in 1:model.T
            flow_val = flows[a, r, t]
            violation_amount = abs(flow_val) - line_limit_fixed
            ub = (flow_val >= 0)
            key = (a, r, t, ub)

            # Update statistics (tracked constraints get tightness recorded every iter)
            if haskey(model.constraint_stats, key)
                stat = model.constraint_stats[key]
                push!(stat.tightness_history, (niter, violation_amount))
                if violation_amount > model.ptdf_tol
                    stat.total_appearances += 1
                    stat.max_violation = max(stat.max_violation, violation_amount)
                end
            elseif violation_amount > model.ptdf_tol
                model.constraint_stats[key] = ConstraintStat(
                    niter, 1, violation_amount, [(niter, violation_amount)]
                )
            end

            # Only add to violation list if not already constrained
            if get(model.tracked_constraints, key, false)
                continue
            end
            if violation_amount > model.ptdf_tol
                push!(violations, (a, r, t, ub, violation_amount))
            end
        end
    end

    return violations
end

# Remove tracked constraints that have been slack for the last stale_k iterations
function prune_stale_constraints!(model::PTDFModel, stale_k::Int)
    to_prune = Tuple{Int,Int,Int,Bool}[]

    for (key, _) in model.tracked_constraints
        stat = get(model.constraint_stats, key, nothing)
        isnothing(stat) && continue
        history = stat.tightness_history
        length(history) < stale_k && continue
        if stat.total_appearances == 1 && all(v <= model.ptdf_tol for (_, v) in history[end-stale_k+1:end])
            push!(to_prune, key)
        end
    end

    for key in to_prune
        (a, r, t, ub) = key
        ckey = ub ? "$(a)_$(r)_$(t)_ub" : "$(a)_$(r)_$(t)_lb"
        if haskey(model.jump_model[:ptdf_flow], ckey)
            delete(model.jump_model, model.jump_model[:ptdf_flow][ckey])
            delete!(model.jump_model[:ptdf_flow], ckey)
        end
        delete!(model.tracked_constraints, key)
    end

    return length(to_prune)
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
        violations = find_violations(model, flows, gamma_values, niter + 1)
        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        
        # Add constraints for top violations
        n_added = add_constraints!(model, sorted_violations)

        stale_k = get(model.data["param"], "stale_k", nothing)
        n_pruned = (isnothing(stale_k) || n_violated == 0) ? 0 : prune_stale_constraints!(model, stale_k)

        if n_added > 0 || n_pruned > 0
            save_tracked_constraints(model.simdir, model, n_violated)
        end

        solved = (n_violated == 0)
        niter += 1
        println("Iteration $niter: Found $n_violated violations, added $n_added constraints, pruned $n_pruned stale")
    end

    # Record solve time
    model.solve_time = time() - t0
    save_solve_time(model.simdir, model.solve_time)
    save_power_injections(model.simdir, model.jump_model, model.data)
    save_constraint_stats(model.simdir, model)
    
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

function save_constraint_stats(simdir, model::PTDFModel)
    output_dir = joinpath(simdir, "output")
    mkpath(output_dir)

    stats = model.constraint_stats
    if isempty(stats)
        return
    end

    # Summary: one row per constraint
    summary_df = DataFrame(
        arc             = [k[1] for k in keys(stats)],
        rep             = [k[2] for k in keys(stats)],
        time            = [k[3] for k in keys(stats)],
        ub              = [k[4] for k in keys(stats)],
        first_iter      = [v.first_iter        for v in values(stats)],
        total_appearances = [v.total_appearances for v in values(stats)],
        max_violation   = [v.max_violation     for v in values(stats)],
        num_iters_tracked = [length(v.tightness_history) for v in values(stats)],
    )
    CSV.write(joinpath(output_dir, "constraint_stats_summary.csv"), summary_df)

    # Tightness history: one row per (constraint, iteration)
    arcs  = Int[]
    reps  = Int[]
    times = Int[]
    ubs   = Bool[]
    iters = Int[]
    vals  = Float64[]
    for (k, v) in stats
        for (iter, val) in v.tightness_history
            push!(arcs,  k[1]); push!(reps,  k[2])
            push!(times, k[3]); push!(ubs,   k[4])
            push!(iters, iter); push!(vals,  val)
        end
    end
    history_df = DataFrame(arc=arcs, rep=reps, time=times, ub=ubs,
                           iteration=iters, violation_amount=vals)
    CSV.write(joinpath(output_dir, "constraint_tightness_history.csv"), history_df)
end

function set_previous_investments!(
    model::PTDFModel,
    gamma_min::Vector,
    s_energy_min::Vector;
    export_csv::Bool=true
)
    """
    Set previous upgrade levels for investments in the PTDFModel.
    This ensures new investments are at least as large as previous investments.
    
    Parameters:
    -----------
    model : PTDFModel
        The model to set previous upgrades in
    gamma_min : Vector
        Previous line investment levels for all lines
    s_energy_min : Vector
        Previous storage investment levels for all nodes
    export_csv : Bool
        Whether to export previous upgrades to CSV (default: true)
    """
    
    E = model.E
    N = model.N
    
    # Validate dimensions
    if length(gamma_min) != E
        throw(ArgumentError("gamma_min length ($(length(gamma_min))) must match number of lines ($E)"))
    end
    
    if length(s_energy_min) != N
        throw(ArgumentError("s_energy_min length ($(length(s_energy_min))) must match number of nodes ($N)"))
    end
    
    # Export previous upgrades if requested
    if export_csv
        min_upgrades_dir = joinpath(model.simdir, "previous_investment_dir")
        export_investments_csv(
            model.data, 
            gamma_min,
            s_energy_min,
            output_dir=min_upgrades_dir
        )
        @info "Exported previous upgrades to $min_upgrades_dir"
    end
    
    # Handle previous line investments
    # Remove existing constraint if it exists
    if haskey(model.jump_model.obj_dict, :min_gamma)
        delete(model.jump_model, model.jump_model[:min_gamma])
        unregister(model.jump_model, :min_gamma)
    end
    
    # Set previous line investment levels
    @constraint(
        model.jump_model,
        min_gamma[a=1:E],
        model.jump_model[:gamma][a] >= gamma_min[a]
    )
    
    num_line_constraints = count(x -> x > 0, gamma_min)
    if num_line_constraints > 0
        @info "Set previous line upgrades for $num_line_constraints lines"
    end
    
    # Handle previous storage investments
    # Remove existing constraint if it exists
    if haskey(model.jump_model.obj_dict, :min_energy)
        delete(model.jump_model, model.jump_model[:min_energy])
        unregister(model.jump_model, :min_energy)
    end
    
    # Set previous storage investment levels
    @constraint(
        model.jump_model,
        min_energy[i=1:N],
        model.jump_model[:s_energy][i] >= s_energy_min[i]
    )
    
    num_storage_constraints = count(x -> x > 0, s_energy_min)
    if num_storage_constraints > 0
        @info "Set previous storage upgrades for $num_storage_constraints nodes"
    end
    
    return nothing
end

function fix_investments!(
    model::PTDFModel, 
    gamma_val::Union{Vector,Nothing}=nothing, 
    s_energy_val::Union{Vector,Nothing}=nothing;
    export_csv::Bool=true
)
    """
    Fix investment decisions in the PTDFModel for Benders decomposition.
    This is used when the model serves as a subproblem.
    
    Removes integer restrictions and bounds from fixed variables to enable
    dual variable access (required for Benders cuts).
    
    Parameters:
    -----------
    model : PTDFModel
        The model to fix investments in
    gamma_val : Vector or Nothing
        Line investment values. If nothing, line investments won't be fixed
    s_energy_val : Vector or Nothing  
        Storage investment values. If nothing, storage investments won't be fixed
    export_csv : Bool
        Whether to export investments to CSV (default: true)
    
    Examples:
    ---------
    # Fix both line and storage investments
    fix_investments!(model, gamma_val, s_energy_val)
    
    # Fix only line investments
    fix_investments!(model, gamma_val, nothing)
    
    # Fix only storage investments
    fix_investments!(model, nothing, s_energy_val)
    """
    
    # Validate that at least one investment type is provided
    if isnothing(gamma_val) && isnothing(s_energy_val)
        @warn "No investment values provided. No constraints will be fixed."
        return
    end
    
    E = model.E
    N = model.N
    
    # Validate dimensions if values are provided
    if !isnothing(gamma_val) && length(gamma_val) != E
        throw(ArgumentError("gamma_val length ($(length(gamma_val))) must match number of lines ($E)"))
    end
    if !isnothing(s_energy_val) && length(s_energy_val) != N
        throw(ArgumentError("s_energy_val length ($(length(s_energy_val))) must match number of nodes ($N)"))
    end
    
    # Export investments if requested
    if export_csv
        current_inv_dir = joinpath(model.simdir, "current_investments")
        export_investments_csv(
            model.data, 
            isnothing(gamma_val) ? zeros(E) : gamma_val,
            isnothing(s_energy_val) ? zeros(N) : s_energy_val,
            output_dir=current_inv_dir
        )
    end
    
    # Handle line investments (gamma)
    if !isnothing(gamma_val)
        gamma = model.jump_model[:gamma]
        
        # Make variables free for dual access
        for a in 1:E
            # Remove integer restriction if present
            if is_integer(gamma[a])
                unset_integer(gamma[a])
            end
            
            # Remove bounds
            if has_lower_bound(gamma[a])
                delete_lower_bound(gamma[a])
            end
            if has_upper_bound(gamma[a])
                delete_upper_bound(gamma[a])
            end
        end
        
        # Remove existing constraint if it exists
        if haskey(model.jump_model.obj_dict, :master_gamma)
            delete(model.jump_model, model.jump_model[:master_gamma])
            unregister(model.jump_model, :master_gamma)
        end
        
        # Fix line investment decisions
        @constraint(
            model.jump_model, 
            master_gamma[a=1:E], 
            gamma[a] == gamma_val[a]
        )
    end
    
    # Handle storage investments (s_energy)
    if !isnothing(s_energy_val)
        s_energy = model.jump_model[:s_energy]
        
        # Make variables free for dual access
        for i in 1:N
            # Remove bounds
            if has_lower_bound(s_energy[i])
                delete_lower_bound(s_energy[i])
            end
            if has_upper_bound(s_energy[i])
                delete_upper_bound(s_energy[i])
            end
        end
        
        # Remove existing constraint if it exists
        if haskey(model.jump_model.obj_dict, :master_energy)
            delete(model.jump_model, model.jump_model[:master_energy])
            unregister(model.jump_model, :master_energy)
        end
        
        # Fix storage investment decisions
        @constraint(
            model.jump_model, 
            master_energy[i=1:N], 
            s_energy[i] == s_energy_val[i]
        )
    end
    
    return nothing
end

function extract_investment_duals(model::PTDFModel)
    """
    Extract duals from investment-fixing constraints for Benders cuts.
    """
    if termination_status(model.jump_model) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        error("Model not solved optimally: $(termination_status(model.jump_model))")
    end
    
    dual_gamma = [dual(model.jump_model[:master_gamma][a]) for a=1:model.E]
    dual_energy = [dual(model.jump_model[:master_energy][i]) for i=1:model.N]
    
    return (dual_gamma, dual_energy)
end

function set_objective!(model::PTDFModel, objective_type::Symbol)
    """
    Set the objective function for the PTDFModel. Assumes that fix_investments! has been called
    objective_type can be :operational or :feasibility
    """
    data = model.data
    R = model.R
    N = model.N
    T = model.T
    G = model.G
    
    operational_weight = get(data["param"], "operational_weight", 1)
    
    if objective_type == :feasibility
        # Minimize unserved energy only
        @objective(model.jump_model, Min,
            sum(
                data["param"]["representative_prob"][r] *
                sum(sum(model.jump_model[:ue][r, i, t] for i=1:N) for t=1:T)
            for r=1:R)
        )
    elseif objective_type == :operational
        # Full operational cost
        @objective(model.jump_model, Min,
            sum(
                data["param"]["representative_prob"][r] *
                (
                    sum(
                        sum(compute_gen_cost(model.jump_model[:pg][r, g, t], data["gen"]["$g"]) for g=1:G) +
                        sum(data["param"]["under_served_penalty"] * model.jump_model[:ue][r, i, t] for i=1:N) + 
                        sum(get(data["param"], "storage_operation_cost", 0.0) * 
                            (model.jump_model[:ch][r,i,t] + model.jump_model[:dis][r,i,t]) for i=1:N)
                    for t=1:T)
                )
            for r=1:R) * operational_weight
        )
    else
        error("Unknown objective type: $objective_type. Use :operational or :feasibility")
    end
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