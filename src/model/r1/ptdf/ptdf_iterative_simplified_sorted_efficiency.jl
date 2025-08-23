using JuMP
using CSV, DataFrames
include("../../../helpers/compute_gen_cost.jl")
include("../storage_candidates/naive_candidates.jl")
include("../rate_a_zero.jl")
include("base_ptdf.jl")
include("ptdf_save_data.jl")

function create_model_r1_ptdf_iterative_simplified_sorted_efficiency(simdir, data::Dict{String, Any}, optimizer; 
    max_ptdf_iterations::Int=256,
    max_ptdf_per_iteration::Int=32,
    ptdf_tol::Float64=1e-6)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Initialize model and basic components
    model = create_base_model(data, optimizer)
    # Set up model extensions for metadata
    model.ext[:solve_metadata] = Dict(
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
        :ptdf_tol => ptdf_tol
    )
    logfile = joinpath(simdir, "gurobi_logfile.log")
    set_optimizer_attribute(model, "LogFile", logfile)
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])

    # Store PTDF matrix in model extension
    model.ext[:PTDF] = do_all_ptdf(data)

    if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
        ptdf_matrix = model.ext[:PTDF]
        ptdf_sparse = map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"]  # retain only significant entries
        ptdf_trimmed = ptdf_matrix .* ptdf_sparse  # zero out small ones
        model.ext[:PTDF] = ptdf_trimmed
    end
    
    # Track which branches already have constraints
    model.ext[:tracked_constraints] = Dict{Tuple{Int,Int,Int,Bool}, Bool}()

    # Initialize rate_a lookup
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    model.ext[:rate_a_nonzero] = Set(parse(Int, x) for x in rate_a_nonzero)
    
    # Add container for flow constraints
    model[:ptdf_flow] = Dict{String, ConstraintRef}()

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    # Optimized warm-start for PTDF with matrix operations
    if isfile(joinpath(simdir, "tracked_constraints.csv"))
        tracked_df = CSV.read(joinpath(simdir, "tracked_constraints.csv"), DataFrame)
        println("Beginning optimized warm-start process")
        
        # Group the constraints by arc once
        arc_groups = groupby(tracked_df, :arc)
        total_constraints = 0
        
        for group in arc_groups
            a = first(group).arc
            line_limit_base = data["branch"]["$a"]["rate_a"]
            cap_increment = get_capacity_increment(data, a)
            ptdf_row = model.ext[:PTDF][a,:]
            
            # Process constraints for this arc in batches
            constraints_batch = collect(eachrow(group))
            batch_size = length(constraints_batch)
            total_constraints += batch_size
            
            # Extract all the (r,t) pairs for this arc
            rep_time_pairs = [(cons.rep, cons.time, cons.ub) for cons in constraints_batch]
            
            # Pre-compute the line limits for all constraints in this group
            line_limit = line_limit_base + model[:gamma][a] * cap_increment
            
            # Create a matrix expression for faster constraint generation
            for (idx, (r, t, ub)) in enumerate(rep_time_pairs)
                # Build the net injection vector for this (r,t) pair
                # We can create this once and reuse it for both upper and lower bounds
                # This is the vectorized equivalent of inner sum:
                # sum(row[i] * (net_injection at bus i) for i in 1:N)
                
                # Define the constraint expression using matrix operations
                flow_expr = sum(ptdf_row[i] * (
                    sum(model[:pg][r,g,t] for g in 1:G if gen_bus_map[g] == i; init=0.0) -
                    data["bus"]["$i"]["load"]["$r"][t] +
                    model[:ue][r,i,t] -
                    model[:ch][r,i,t] +
                    model[:dis][r,i,t]
                ) for i in 1:N)
                
                # Add constraint that should be added depending on ub (upper bound)
                if ub
                    model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(model, flow_expr <= line_limit)
                else
                    model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(model, flow_expr >= -line_limit)
                end
                
                # Track the constraint
                model.ext[:tracked_constraints][(a,r,t,ub)] = true
            end
            
            # println("Progress: Cumulatively added $total_constraints constraints")
        end
        
        println("Warm-started with $total_constraints constraints")
    end

    # Begin lazy PTDF loop
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < model.ext[:solve_metadata][:max_ptdf_iterations]
        # Solve model
        optimize!(model)
        
        # Exit if not solved optimally
        st = termination_status(model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(model[:pg])
        ue_values = value.(model[:ue])
        ch_values = value.(model[:ch])
        dis_values = value.(model[:dis])
        gamma_values = value.(model[:gamma])
        
        # Compute all flows at once
        flows = zeros(length(model.ext[:rate_a_nonzero]), R, T)
        compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, data, model.ext[:PTDF])

        # Add violations prioritizing by size
        violations = []
        
        # Keep track of all tracked combinations for this branch to avoid lookups
        for a in model.ext[:rate_a_nonzero]
            
            # Pre-compute for this branch
            ptdf_row = model.ext[:PTDF][a,:]
            cap_increment = get_capacity_increment(data, a)
            line_limit_base = data["branch"]["$a"]["rate_a"]
            line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment
            
            for r in 1:R, t in 1:T
                flow_val = flows[a, r, t]
                violation_amount = abs(flow_val) - line_limit_fixed
                ub = (flow_val >= 0)

                if get(model.ext[:tracked_constraints], (a, r, t, ub), false)
                    continue
                end
        
                if violation_amount > model.ext[:solve_metadata][:ptdf_tol]
                    push!(violations, (a, r, t, ub, violation_amount))
                end
            end
        end

        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        n_added = 0

        # Add up to the limit
        max_to_add = model.ext[:solve_metadata][:max_ptdf_per_iteration]
        for (a, r, t, ub, _) in Iterators.take(sorted_violations, max_to_add)

            ptdf_row = model.ext[:PTDF][a,:]
            cap_increment = get_capacity_increment(data, a)
            line_limit_base = data["branch"]["$a"]["rate_a"]
            line_limit_expr = line_limit_base + model[:gamma][a] * cap_increment
        
            # Define the flow expression only once
            flow_expr = sum(ptdf_row[i] * (
                sum(model[:pg][r, g, t] for g in 1:G if gen_bus_map[g] == i; init=0.0)
                - data["bus"]["$i"]["load"]["$r"][t]
                + model[:ue][r, i, t] - model[:ch][r, i, t] + model[:dis][r, i, t]
            ) for i in 1:N)
        
            # Add both constraints
            if ub
                model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(model, flow_expr <= line_limit_expr)
            else
                model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(model, flow_expr >= -line_limit_expr)
            end
        
            model.ext[:tracked_constraints][(a, r, t, ub)] = true
            n_added += 1
        end

        if n_added > 0
            save_tracked_constraints(simdir, model, n_violated)
        end
        
        solved = (n_violated == 0)
        niter += 1
        # println("approx. progress $(checked_count / E)")
        println("Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    solve_time = time() - t0
    model.ext[:solve_metadata][:solve_time] = solve_time

    save_solve_time(simdir, solve_time)
    save_power_injections(simdir, model, data)
    return model
end

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
    
    # generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    JuMP.@variable(model, pg[r in 1:R, g in 1:G, t in 1:T])

    for g in 1:G
        gen = data["gen"]["$g"]
        is_renewable = gen["gen_type"] ∈ data["param"]["renewable_types"]
        is_foreign = gen["gen_type"]  ∈ ["foreign"]
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
        JuMP.@variable(model, 0 <= gamma[a=1:E] <= K)
    else
        # Use integer variables for standard model
        JuMP.@variable(model, 0 <= gamma[a=1:E] <= K, Int) # investment level of capacity upgrade
    end
    
    JuMP.@variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0) # under-served energy at bus
    # JuMP.@variable(model, pf[r=1:R, a=1:E, t=1:T]) # branch flows
    # JuMP.@variable(model, s_power[i=1:N] >= 0) # power rating of storage
    JuMP.@variable(model, s_energy[i=1:N] >= 0) # energy rating of storage
    JuMP.@variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0) # state of charge of storage
    JuMP.@variable(model, ch[r=1:R, i=1:N, t=1:T] >= 0) # charging of storage
    JuMP.@variable(model, dis[r=1:R, i=1:N, t=1:T] >= 0) # discharging of storage
    # JuMP.@variable(model, sigma[i=1:N], Bin) # binary variable for installation of storage

    # Add all non-PTDF constraints
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

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # soc over time constraint
    JuMP.@constraint(model, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - dis[r,i,t] / data["param"]["bess_efficiency"]
    )

    # OPTIONAL: soc init and end constraint
    JuMP.@constraint(model,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] + ch[r,i,1] * data["param"]["bess_efficiency"] - dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    JuMP.@constraint(model,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i]
    )

    # soc energy rating constraint
    JuMP.@constraint(model, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i]
    )

    # energy rating only if storage installed
    JuMP.@constraint(model, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= data["param"]["max_energy_rating"]
    )

    # charge/discharge must be constrained by power rating
    JuMP.@constraint(model,
        charge_discharge_lb[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] <= s_energy[i] / 4
    )
    JuMP.@constraint(model,
        charge_discharge_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= s_energy[i] / 4
    )

    # OPTIONAL: precompute maximum nodal injections
    if get(data["param"], "nodecap", false)
        nodal_caps = Dict{Int, Float64}()
        for i in 1:N
            max_outflow = sum(data["branch"]["$a"]["rate_a"] + get_capacity_increment(data, a) * K for a in data["arcs_from"]["$i"])
            
            max_charge = min(data["param"]["max_power_rating"], data["param"]["max_energy_rating"] / 4)
            nodal_caps[i] = max_outflow + max_charge                
        end

        JuMP.@constraint(model,
            nodal_injection_cap[r in 1:R, i in 1:N, t in 1:T],
            sum(model[:pg][r, g, t] for g in 1:G if gen_bus_map[g] == i) <=
            data["bus"]["$i"]["load"]["$r"][t] + nodal_caps[i]
        )
    end
    
    # OPTIONAL: CANDIDATE STORAGE LOCATIONS ONLY
    cand_file = joinpath(simdir, "cand.json")
    if isfile(cand_file)
        cand_dict = JSON.parsefile(cand_file)
        candidates = Set(string.(cand_dict["candidates"]))

        all_busses = Set(keys(data["bus"]))
        non_candidates = setdiff(all_busses, union(candidates, nonzero_storage_nodes))
        non_candidates = Set(parse(Int, x) for x in non_candidates)

        for i in non_candidates
            fix(s_energy[i], 0; force = true)
        end

        println("Number of candidates in Gurobi model: $(length(candidates))")
    end

    # Objective
    operational_weight = get(data["param"], "operational_weight", 1)
    @objective(model, Min,
        sum(s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
        sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g in 1:G) +
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i in 1:N) +
                    sum(get(data["param"], "storage_operation_cost", 0.0) * (ch[r,i,t] + dis[r,i,t]) for i=1:N)
                for t in 1:T)
            )
        for r in 1:R) * operational_weight
    )
    
    return model
end

function compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, data, PTDF)
    # Preallocate net injection vector
    n_buses = length(data["bus"])
    net_injections = zeros(n_buses)  # PTDF matrix is (n_branches × n_buses)
    
    # For each timestep
    for r in 1:size(pg_values, 1), t in 1:size(pg_values, 3)
        fill!(net_injections, 0.0)

        # Compute net injections at each bus
        for i in 1:n_buses 
            # Sum generation at bus i
            net_injections[i] = sum(pg_values[r,g,t] 
                                    for g in 1:length(data["gen"]) 
                                    if data["gen"]["$g"]["gen_bus"] == i; 
                                    init=0.0)
            # Subtract load, add unserved energy, subtract charging, add discharging
            net_injections[i] -= data["bus"]["$i"]["load"]["$r"][t]
            net_injections[i] += ue_values[r,i,t]
            net_injections[i] -= ch_values[r,i,t]
            net_injections[i] += dis_values[r,i,t]
        end
        
        # Compute flows for all branches at once using PTDF
        flows[:,r,t] = PTDF * net_injections
    end
end

function add_prev_upgrades(model, data, gamma, s_energy, prev_dir; tolerance=1e-5)
    E = length(data["branch"])
    N = length(data["bus"])

    trans_file = joinpath(prev_dir, "line_investments.csv")
    storage_file = joinpath(prev_dir, "storage_investments.csv")

    trans_df = CSV.read(trans_file, DataFrame)
    storage_df = CSV.read(storage_file, DataFrame)

    @constraint(model,
        old_gamma[a in 1:E, trans_df[a, :Upgrade_Lvl] > tolerance],
        gamma[a] >= trans_df[a, :Upgrade_Lvl]
    )
    @constraint(model,
        old_s_energy[i in 1:N, storage_df[i, :Storage_Energy] > tolerance],
        s_energy[i] >= storage_df[i, :Storage_Energy]
    )
end