function save_tracked_constraints(simdir, tracked_constraints, n_violated)
    # Convert dictionary to DataFrame
    df = DataFrame(
        arc = [k[1] for k in keys(tracked_constraints)],
        rep = [k[2] for k in keys(tracked_constraints)],
        time = [k[3] for k in keys(tracked_constraints)],
        ub = [k[4] for k in keys(tracked_constraints)],
        tracked = collect(values(tracked_constraints)),
        violations_left = [n_violated for i in keys(tracked_constraints)]
    )

    # Sort DataFrame by arc first, then rep, then time
    # sort!(df, [:arc, :rep, :time])

    # Save to CSV
    CSV.write(joinpath(simdir, "output", "tracked_constraints.csv"), df)
end

function create_subproblem_model(simdir, data, tracked_constraints)
    """
    Creates the subproblem model with all variables and constraints except the objective function.
    Returns the model ready for objective setting and solving.
    """

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Initialize model
    sub_optimizer = Gurobi.Optimizer
    sub = JuMP.Model(sub_optimizer)
    set_optimizer_attribute(sub, "LogFile", joinpath(simdir, "gurobi_sub_logfile.log"))
    set_optimizer_attribute(sub, "MIPGap", data["param"]["mip_gap"])
    set_optimizer_attribute(sub, "Method", 1)      # Force dual simplex
    set_optimizer_attribute(sub, "Crossover", 0)   # Disable crossover

    # Store additional data in model extensions
    sub.ext[:data] = data
    sub.ext[:simdir] = simdir
    
    # Compute PTDF matrix
    sub.ext[:PTDF] = do_all_ptdf(data)

    if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
        ptdf_matrix = sub.ext[:PTDF]
        ptdf_sparse = map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"]  # retain only significant entries
        ptdf_trimmed = ptdf_matrix .* ptdf_sparse  # zero out small ones
        sub.ext[:PTDF] = ptdf_trimmed
    end

    # Initialize tracked constraints with the passed dictionary
    sub.ext[:tracked_constraints] = deepcopy(tracked_constraints)
    
    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    sub.ext[:rate_a_nonzero] = Set(parse(Int, x) for x in rate_a_nonzero)

    # Container for flow constraints
    sub[:ptdf_flow] = Dict{String, ConstraintRef}()

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))
    sub.ext[:gen_bus_map] = gen_bus_map

    #
    #   I. Variables
    #
    
    # generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    JuMP.@variable(sub, pg[r in 1:R, g in 1:G, t in 1:T])

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

    # Under-served energy
    @variable(sub, ue[r=1:R, i=1:N, t=1:T] >= 0)
    
    # Storage operation
    @variable(sub, ch[r=1:R, i=1:N, t=1:T] >= 0)  # Charging/discharging
    @variable(sub, dis[r=1:R, i=1:N, t=1:T] >= 0)
    @variable(sub, soc[r=1:R, i=1:N, t=1:T] >= 0)  # State of charge

    # Fix investment variables based on master problem decisions
    @variable(sub, gamma[a=1:E])
    @variable(sub, s_energy_int[i=1:N])

    #
    #   II. Constraints
    #

    # FIX INVESTMENT DECISIONS FROM MASTER PROBLEM (DUALS USED FOR CUTS)
    # @constraint(sub, master_gamma[a=1:E], gamma[a] == gamma_val[a])
    # @constraint(sub, master_energy[i=1:N], s_energy_int[i] == s_energy_int_val[i])

    # Global power balance
    @constraint(sub, 
        power_balance[r=1:R, t=1:T],
        sum(pg[r,g,t] for g=1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i=1:N)
        + sum(ue[r,i,t] for i=1:N)
        - sum(ch[r,i,t] for i=1:N)
        + sum(dis[r,i,t] for i=1:N)
        == 0
    )

    # Storage constraints

    # SOC evolution
    @constraint(sub, 
        soc_over_time[r=1:R, i=1:N, t=2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - dis[r,i,t] / data["param"]["bess_efficiency"]
    )
    
    # Initial and final SOC
    @constraint(sub,
        soc_start[r=1:R, i=1:N],
        soc[r,i,1] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy_int[i] * data["param"]["storage_energy_size"] + ch[r,i,1] * data["param"]["bess_efficiency"] - dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    @constraint(sub,
        soc_end[r=1:R, i=1:N],
        soc[r,i,T] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy_int[i] * data["param"]["storage_energy_size"]
    )

    # SOC limits
    @constraint(sub, 
        soc_energy_ub[r=1:R, i=1:N, t=1:T],
        soc[r,i,t] <= s_energy_int[i] * data["param"]["storage_energy_size"]
    )
    
    # Charging/discharging limits
    @constraint(sub,
        charge_discharge_lb[r=1:R, i=1:N, t=1:T],
        dis[r,i,t] <= s_energy_int[i] * data["param"]["storage_energy_size"] / 4
    )
    @constraint(sub,
        charge_discharge_ub[r=1:R, i=1:N, t=1:T],
        ch[r,i,t] <= s_energy_int[i] * data["param"]["storage_energy_size"] / 4
    )

    # Apply warm-start using tracked constraints
    if !isempty(tracked_constraints)
        println("[DEBUG] Beginning optimized warm-start for subproblem")
        
        # Convert dictionary to DataFrame for grouping
        arc_list = [k[1] for k in keys(tracked_constraints)]
        rep_list = [k[2] for k in keys(tracked_constraints)]
        time_list = [k[3] for k in keys(tracked_constraints)]
        ub_list = [k[4] for k in keys(tracked_constraints)]
        tracked_df = DataFrame(arc=arc_list, rep=rep_list, time=time_list, ub=ub_list)
        
        # Group by arc for efficient processing
        arc_groups = groupby(tracked_df, :arc)
        total_constraints = 0
        
        for group in arc_groups
            a = first(group).arc
            line_limit_base = data["branch"]["$a"]["rate_a"]
            cap_increment = get_capacity_increment(data, a)
            ptdf_row = sub.ext[:PTDF][a,:]
            
            constraints_batch = collect(eachrow(group))
            batch_size = length(constraints_batch)
            total_constraints += batch_size
            
            # Get all (r,t) pairs for this arc
            rep_time_pairs = [(row.rep, row.time, row.ub) for row in constraints_batch]

            # Pre-compute the line limits for all constraints in this group
            line_limit = line_limit_base + sub[:gamma][a] * cap_increment
            
            # Generate constraints efficiently
            for (r, t, ub) in rep_time_pairs
                # Create flow expression once and reuse
                flow_expr = sum(ptdf_row[i] * (
                    sum(sub[:pg][r,g,t] for g in 1:G if gen_bus_map[g] == i; init=0.0) -
                    data["bus"]["$i"]["load"]["$r"][t] +
                    sub[:ue][r,i,t] -
                    sub[:ch][r,i,t] +
                    sub[:dis][r,i,t]
                ) for i in 1:N)
                
                # Add constraint that should be added depending on ub (upper bound)
                if ub
                    sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit)
                else
                    sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit)
                end
            end
        end
        
        println("[DEBUG] Warm-started subproblem with $total_constraints constraints")
    end

    return sub
end

function set_operational_objective!(sub)
    """
    Sets the operational cost objective function on the subproblem model.
    """
    data = sub.ext[:data]
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    
    operational_weight = get(data["param"], "operational_weight", 1)
    @objective(sub, Min,
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(compute_gen_cost(sub[:pg][r, g, t], data["gen"]["$g"]) for g=1:G) +
                    sum(data["param"]["under_served_penalty"] * sub[:ue][r, i, t] for i=1:N) + 
                    sum(get(data["param"], "storage_operation_cost", 0.0) * (sub[:ch][r,i,t] + sub[:dis][r,i,t]) for i=1:N)
                for t=1:T)
            )
        for r=1:R) * operational_weight
    )
end

function set_underserved_objective!(sub)
    """
    Sets the underserved energy objective function on the subproblem model.
    """
    data = sub.ext[:data]
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    T = data["param"]["num_hours"]
    
    @objective(sub, Min,
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(sub[:ue][r, i, t] for i=1:N)
                for t=1:T)
            )
        for r=1:R)
    )
end

function set_sub_investments!(sub, y_val)
    data = sub.ext[:data]
    N = length(data["bus"])
    E = length(data["branch"])

    gamma = sub[:gamma]
    s_energy_int = sub[:s_energy_int]

    # Check and remove existing master_gamma constraints if they exist
    if haskey(sub.obj_dict, :master_gamma)
        delete(sub, sub[:master_gamma])
        unregister(sub, :master_gamma)
    end
    
    # Check and remove existing master_energy constraints if they exist
    if haskey(sub.obj_dict, :master_energy)
        delete(sub, sub[:master_energy])
        unregister(sub, :master_energy)
    end

    # FIX INVESTMENT DECISIONS FROM MASTER PROBLEM (DUALS USED FOR CUTS)
    gamma_val, s_energy_int_val = y_val
    @constraint(sub, master_gamma[a=1:E], gamma[a] == gamma_val[a])
    @constraint(sub, master_energy[i=1:N], s_energy_int[i] == s_energy_int_val[i])
end

function solve_subproblem_ptdf!(sub; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6, logging=nothing)
    """
    Solves the subproblem using lazy PTDF constraints.
    The model should already be created and have an objective function set.
    """
    
    if !isnothing(logging)
        println("[DEBUG] $logging: Beginning subproblem solve")
    end
    
    # Get data from model extensions
    data = sub.ext[:data]
    simdir = sub.ext[:simdir]
    gen_bus_map = sub.ext[:gen_bus_map]
    
    # Set up the PTDF metadata
    sub.ext[:solve_metadata] = Dict(
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
        :ptdf_tol => ptdf_tol
    )
    
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])

    # Solve using lazy PTDF constraints
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < max_ptdf_iterations
        # Solve model
        println("[DEBUG] Solving PTDF subproblem for constraint discovery: subiteration $niter")
        optimize!(sub)
        
        # Exit if not solved optimally
        st = termination_status(sub)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(sub[:pg])
        ue_values = value.(sub[:ue])
        ch_values = value.(sub[:ch])
        dis_values = value.(sub[:dis])
        gamma_values = value.(sub[:gamma])

        # Compute all flows at once
        flows = zeros(length(sub.ext[:rate_a_nonzero]), R, T)
        compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, data, sub.ext[:PTDF])

        # Add violations prioritizing by size
        violations = []

        # Keep track of all tracked combinations for this branch to avoid lookups
        for a in sub.ext[:rate_a_nonzero]
            
            # Pre-compute for this branch
            ptdf_row = sub.ext[:PTDF][a,:]
            cap_increment = get_capacity_increment(data, a)
            line_limit_base = data["branch"]["$a"]["rate_a"]
            line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment
            
            for r in 1:R, t in 1:T
                flow_val = flows[a, r, t]
                violation_amount = abs(flow_val) - line_limit_fixed
                ub = (flow_val >= 0)

                if get(sub.ext[:tracked_constraints], (a, r, t, ub), false)
                    continue
                end
        
                if violation_amount > ptdf_tol
                    push!(violations, (a, r, t, ub, violation_amount))
                end
            end
        end

        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        n_added = 0

        # Add up to the limit
        for (a, r, t, ub, _) in Iterators.take(sorted_violations, max_ptdf_per_iteration)

            ptdf_row = sub.ext[:PTDF][a,:]
            cap_increment = get_capacity_increment(data, a)
            line_limit_base = data["branch"]["$a"]["rate_a"]
            line_limit_expr = line_limit_base + sub[:gamma][a] * cap_increment
        
            # Define the flow expression only once
            flow_expr = sum(ptdf_row[i] * (
                sum(sub[:pg][r, g, t] for g in 1:G if gen_bus_map[g] == i; init=0.0)
                - data["bus"]["$i"]["load"]["$r"][t]
                + sub[:ue][r, i, t] - sub[:ch][r, i, t] + sub[:dis][r, i, t]
            ) for i in 1:N)
        
            # Add both constraints
            if ub
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit_expr)
            else
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit_expr)
            end
        
            sub.ext[:tracked_constraints][(a, r, t, ub)] = true
            n_added += 1
        end

        if n_added > 0
            save_tracked_constraints(simdir, sub, n_violated)
        end
        solved = (n_violated == 0)
        niter += 1

        println("[DEBUG] PTDF Subproblem Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    solve_time = time() - t0
    sub.ext[:solve_metadata][:solve_time] = solve_time

    return sub
end

function extract_subproblem_duals(sub)
    """
    Extracts duals from the solved subproblem for Benders cuts.
    """
    # Extract duals for Benders cuts
    if termination_status(sub) ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        E = length(sub.ext[:data]["branch"])
        N = length(sub.ext[:data]["bus"])
        
        dual_gamma = [dual(sub[:master_gamma][a]) for a=1:E]
        dual_energy = [dual(sub[:master_energy][i]) for i=1:N]
        
        return (dual_gamma, dual_energy)
    else
        error("Subproblem failed to solve optimally: $(termination_status(sub))")
    end
end