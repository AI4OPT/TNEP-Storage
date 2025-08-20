using JuMP
using Gurobi
using Serialization

function define_master_ptdf(superdir, master_data::Dict{String, Any}, date_weights::Dict{Int, Tuple{String, Float64}})
    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)
    
    if !get(master_data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", master_data["param"]["mip_gap"])
    end

    # Initialize sets
    RR = master_data["param"]["num_representatives"]
    N = length(master_data["bus"])
    E = length(master_data["branch"])
    T = master_data["param"]["num_hours"]
    G = length(master_data["gen"])
    K = master_data["param"]["num_cap_upgrades_max"]

    # Create extension for tracking iterations of y_raw, y_eval, and y_core
    master.ext[:y_raw] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_eval] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_core] = Vector{Vector{Vector{Float64}}}()

    # Create extensions for dynamic lambda adjustments
    master.ext[:gap] = Vector{Float64}()
    master.ext[:stabilization_lambda] = Vector{Float64}()
    master.ext[:stabilization_lambda_decr] = Vector{Float64}()
    master.ext[:iter] = 0

    # Create extension for date_weights
    master.ext[:date_weights] = date_weights

    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(master_data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => master_data["gen"]["$g"]["gen_bus"] for g in keys(master_data["gen"]))

    #
    #   I. Variables
    #

    # investment level of capacity upgrade
    if haskey(master_data["param"], "relaxed_first_stage") && master_data["param"]["relaxed_first_stage"] == true
        JuMP.@variable(master, 0 <= gamma[a=1:E] <= K)
    else
        JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)
    end

    # binary variable for installation of storage
    # JuMP.@variable(master, sigma[i=1:N], Bin)

    # power rating of storage
    JuMP.@variable(master, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(master, s_energy[i=1:N] >= 0)

    # subproblem objective(s)
    JuMP.@variable(master, theta[r=1:RR] >= 0)

    #
    #   II. Constraints
    #
    candidate_branches = get_candidate_branches(master_data)
    
    # Check if previous investments exist
    if haskey(master_data["param"], "previous_investment_dir")
        prev_dir = master_data["param"]["previous_investment_dir"]
        add_prev_upgrades(master, master_data, gamma, s_energy, prev_dir, candidate_branches)
    end

    # Separately, force non-candidates to zero
    if candidate_branches !== nothing
        @constraint(master, 
            force_zero[a in rate_a_nonzero; !(a in candidate_branches)],
            gamma[a] == 0
        )
    end

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(master, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # max energy rating
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= master_data["param"]["max_energy_rating"]
    )

    # max power rating
    JuMP.@constraint(master, 
        installed_power_ub[i in 1:N],
        s_power[i] <= master_data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    JuMP.@constraint(master, 
        short_duration[i in 1:N],
        s_energy[i] == 4.0 * s_power[i]
    )

    embed_date = master_data["param"]["embed_date"]
    simdir = joinpath(superdir, embed_date)

    data, embed_obj_expr = embed_subproblem(master, simdir)
    initialize_master_ptdf_tracked(master, simdir, data)

    #
    #   III. Objective
    #
    obj_expr = (sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            sum(theta[r] * date_weights[r][2] for r in 1:RR))

    embed_date_index = findfirst(r -> date_weights[r][1] == embed_date, 1:RR)
    JuMP.@constraint(master, theta[embed_date_index] == embed_obj_expr)

    JuMP.@objective(master, Min, obj_expr)
    
    y = gamma, s_power, s_energy
    return master, y, theta
end

function solve_master!(model, simdir, data)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    
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
end

function initialize_master_ptdf_tracked(model, simdir, data)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    # Set up model extensions for metadata
    model.ext[:solve_metadata] = Dict(
        :max_ptdf_iterations => 256,
        :max_ptdf_per_iteration => 32,
        :ptdf_tol => 1e-6
    )
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

    # Add container for flow constraints
    model.ext[:rate_a_nonzero] = Set(parse(Int, x) for x in rate_a_nonzero)
    model[:ptdf_flow] = Dict{String, ConstraintRef}()

    # Optimized warm-start for PTDF with matrix operations
    filename = joinpath(simdir, "output", "tracked_constraints.csv")
    if isfile(filename) || isfile(joinpath(simdir, "tracked_constraints.csv"))
        if isfile(filename)
            tracked_df = CSV.read(filename, DataFrame)
        else
            tracked_df = CSV.read(joinpath(simdir, "tracked_constraints.csv"), DataFrame)
        end
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
end

function embed_subproblem(model, simdir)
    data = JSON.parsefile(joinpath(simdir, "data.json"))

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    #
    #   I. Variables
    #

    s_energy = model[:s_energy]
    s_power = model[:s_power]
    
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

    JuMP.@variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0) # under-served energy at bus
    JuMP.@variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0) # state of charge of storage
    JuMP.@variable(model, ch[r=1:R, i=1:N, t=1:T] >= 0) # charging of storage
    JuMP.@variable(model, dis[r=1:R, i=1:N, t=1:T] >= 0) # discharging of storage

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

    # charge/discharge must be constrained by power rating
    JuMP.@constraint(model,
        charge_discharge_lb[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] <= s_power[i]
    )
    JuMP.@constraint(model,
        charge_discharge_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= s_power[i]
    )

    operational_weight = get(data["param"], "operational_weight", 1)
    embed_obj_expr = sum(sum(
                    sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g in 1:G) +
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i in 1:N) +
                    sum(get(data["param"], "storage_operation_cost", 0.0) * (ch[r,i,t] + dis[r,i,t]) for i=1:N)
                for t in 1:T)
        for r in 1:R) * operational_weight

    return data, embed_obj_expr
end

function update_master_objective!(superdir, master, data, y, theta)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    date_weights = master.ext[:date_weights]

    #
    #   III. Objective
    #
    gamma, s_power, s_energy = y

    # Start with the base objective terms
    obj_expr = (sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            sum(theta[r] * date_weights[r][2] for r in 1:R))

    # L2 (quadratic) regularization
    if haskey(data["param"], "reg_penalty")
        gamma_reg, s_power_reg, s_energy_reg = compute_regularization_point(superdir, master, data)
        obj_expr += data["param"]["reg_penalty"] * sum((s_power[i] - s_power_reg[i])^2 for i in 1:N)
        if haskey(data["param"], "trans_reg_penalty")
            obj_expr += data["param"]["trans_reg_penalty"] * sum((gamma[a] - gamma_reg[a])^2 for a in 1:E)
        end
    end

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)
end

function add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
    # Unpack y variables and duals
    gamma, s_power, s_energy = y
    dual_gamma, dual_power = duals
    gamma_val, s_power_val, s_energy_val = y_val

    # Add Benders cut
    @constraint(master, 
        theta >= phi_val + 
        sum(dual_gamma[a] * (gamma[a] - gamma_val[a]) for a=1:length(gamma)) +
        sum(dual_power[i] * (s_power[i] - s_power_val[i]) for i=1:length(s_power))
    )
end

function add_prev_upgrades(master, data, gamma, s_energy, prev_dir, candidate_branches=nothing; tolerance=1e-5)
    E = length(data["branch"])
    N = length(data["bus"])

    trans_file = joinpath(prev_dir, "line_investments.csv")
    storage_file = joinpath(prev_dir, "storage_investments.csv")
    trans_df = CSV.read(trans_file, DataFrame)
    storage_df = CSV.read(storage_file, DataFrame)

    if candidate_branches !== nothing
        valid_branches = [a for a in 1:E if (a in candidate_branches) && (trans_df[a, :Upgrade_Lvl] > tolerance)]
    else
        # Original behavior if no candidate set provided
        valid_branches = [a for a in 1:E if trans_df[a, :Upgrade_Lvl] > tolerance]
    end

    @constraint(master,
        old_gamma[a in valid_branches],
        gamma[a] >= trans_df[a, :Upgrade_Lvl]
    )
    @constraint(master,
        old_s_energy[i in 1:N, storage_df[i, :Storage_Energy] > tolerance],
        s_energy[i] >= storage_df[i, :Storage_Energy]
    )
end

function get_candidate_branches(data, tolerance=1e-5)
    if !get(data["param"], "line_candidates", false)
        return nothing
    end
    
    dates = data["param"]["dates"] # Read data params
    initial_optima_dir = data["param"]["initial_optima_dir"]
    
    branch_rep_upgrades = Dict{Int, Vector{Float64}}() # Initialize storage for all branch data
    
    for branch_idx in 1:length(data["branch"]) # Initialize with empty vectors for each branch
        branch_rep_upgrades[branch_idx] = Float64[]
    end
    
    for a_date in dates # Collect rep day data
        rep_file = joinpath(initial_optima_dir, a_date, "output", "line_investments.csv")
        rep = CSV.read(rep_file, DataFrame)
        
        for row in eachrow(rep) # Add each branch's upgrade level to its collection
            push!(branch_rep_upgrades[row.Branch_Index], row.Upgrade_Lvl)
        end
    end

    filtered_dict = Dict(k => v for (k, v) in branch_rep_upgrades if any(abs(val) > tolerance for val in v))
    return collect(keys(filtered_dict))
end



