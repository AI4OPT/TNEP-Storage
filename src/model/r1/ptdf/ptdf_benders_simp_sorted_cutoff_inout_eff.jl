using JuMP
using Gurobi
using Dualization
using CSV, DataFrames
using DataStructures
include("../run_model.jl")
include("../../../helpers/compute_gen_cost.jl")
include("../rate_a_zero.jl")
include("base_ptdf.jl")
include("ptdf_save_data.jl")
# include("ptdf_benders_transport_flow_subproblem.jl")
include("ptdf_core_point.jl")
include("ptdf_benders_pareto.jl")
include("ptdf_benders_save_load_master.jl")

function define_master_ptdf(data::Dict{String, Any})
    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)
    
    if !get(data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])
    end

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Create extension for tracking iterations of y_trust
    master.ext[:y_trust] = Vector{Vector{Vector{Float64}}}()
    master.ext[:trust_transmission] = Vector{Float64}()
    master.ext[:trust_sizing] = Vector{Float64}()
    master.ext[:trust_siting] = Vector{Float64}()

    master.ext[:gap] = Vector{Float64}()
    master.ext[:stabilization_lambda] = Vector{Float64}()
    master.ext[:stabilization_lambda_decr] = Vector{Float64}()
    master.ext[:iter] = 0
    master.ext[:data] = data
    master.ext[:total_ue] = [Inf]
    master.ext[:total_obj] = [Inf]

    # Get incidence matrix from extension
    incidence_matrix = master.ext[:INCIDENCE] = do_all_incidence(data)

    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))
    
    # Branch metadata
    from_bus = Dict(parse(Int, a) => data["branch"]["$a"]["f_bus"] for a in keys(data["branch"]))
    to_bus = Dict(parse(Int, a) => data["branch"]["$a"]["t_bus"] for a in keys(data["branch"]))

    #
    #   I. Variables
    #

    # investment level of capacity upgrade
    JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)

    # binary variable for installation of storage
    # JuMP.@variable(master, sigma[i=1:N], Bin)

    # energy rating of storage
    JuMP.@variable(master, s_energy_int[i=1:N] >= 0, Int)

    # subproblem objective(s)
    JuMP.@variable(master, theta >= 0)

    # subproblem load shed
    JuMP.@variable(master, ue_sum >= 0)

    master.ext[:warmstart] = get(data["param"], "warmstart", false)
    master.ext[:stabilization] = get(data["param"], "stabilization", false)

    if master.ext[:stabilization] == "trust_region"
        JuMP.@variable(master, abs_diff[1:N] >= 0)
        JuMP.@variable(master, trans_abs_diff[1:E] >= 0)
    end

    #
    #   II. Constraints
    #    
    """
    JuMP.@constraint(master, 
        no_load_shed,
        ue_sum <= 0
    )
    """

    """
    # Check if previous investments exist
    if haskey(data["param"], "previous_investment_dir")
        prev_dir = data["param"]["previous_investment_dir"]
        add_prev_upgrades(master, data, gamma, s_energy, prev_dir, nothing)
    end
    """

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(master, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # energy rating maximum
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy_int[i] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
    )

    # storage candidates
    if get(data["param"], "storage_cand", false)
        cand_file = joinpath(data["param"]["core_point_simdir"], "output", "storage_investments.csv")
        storage_df = CSV.read(cand_file, DataFrame)
        storage_upgrades = storage_df[:, :Storage_Energy]
        storage_non_indices = findall(x -> x == 0, storage_upgrades)

        JuMP.@constraint(master, 
            energy_cand[i in storage_non_indices],
            s_energy_int[i] == 0
        )
    end

    # line candidates
    if get(data["param"], "line_cand", false)
        cand_file = joinpath(data["param"]["core_point_simdir"], "output", "line_investments.csv")
        line_df = CSV.read(cand_file, DataFrame)
        line_upgrades = line_df[:, :Upgrade_Lvl]
        line_non_indices = findall(x -> x == 0, line_upgrades)

        JuMP.@constraint(master, 
            line_cand[a in line_non_indices],
            gamma[a] == 0
        )
    end

    #
    #   III. Objective
    #
    # Start with the base objective terms
    obj_expr = sum(s_energy_int[i] * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            theta + ue_sum * 365 * data["param"]["under_served_penalty"]
    y = gamma, s_energy_int

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)

    return master, y, theta
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

function add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
    # Unpack y variables and duals
    gamma, s_energy_int = y
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_int_val = y_val

    # Add Benders cut
    @constraint(master, 
        theta >= phi_val + 
        sum(dual_gamma[a] * (gamma[a] - gamma_val[a]) for a=1:length(gamma)) +
        sum(dual_energy[i] * (s_energy_int[i] - s_energy_int_val[i]) for i=1:length(s_energy_int))
    )
end

function benders_iteration_ptdf(simdir, master, y, theta, data, max_iterations=100000, tolerance=0.01)
    max_iterations=100000
    tolerance=0.01
    converged = false
    iter = master.ext[:iter]
    ue_sum = master[:ue_sum]

    # Initialize tracked constraints dictionary to store across iterations
    tracked_constraints = Dict{Tuple{Int,Int,Int,Bool}, Bool}()

    # Load previous constraints if available
    if isfile(joinpath(simdir, "tracked_constraints.csv"))
        tracked_df = CSV.read(joinpath(simdir, "tracked_constraints.csv"), DataFrame)
        println("[DEBUG] Loading $(nrow(tracked_df)) previously tracked constraints")
        
        for row in eachrow(tracked_df)
            if row.tracked
                tracked_constraints[(row.arc, row.rep, row.time, row.ub)] = true
            end
        end
    end

    # initialize subproblem
    sub = create_subproblem_model(simdir, data, tracked_constraints)

    while !converged && iter < max_iterations
        # Solve master problem
        println("[DEBUG] Solving master problem: iteration $iter")
        add_box_step!(master, 1)
        add_trust_region!(master)
        optimize!(master)
        gamma_val = value.(master[:gamma])
        s_energy_int_val = value.(master[:s_energy_int])
        theta_val = value.(master[:theta])

        # Save lower bound (no trust region)
        master_obj = objective_value(master) 

        line_file = "sim/r1/PSCC2026/nobenders/reformulated/discrete/2035pc_om_d/2016-08-11/output/line_investments.csv"
        # line_inv = CSV.read(line_file, DataFrame)
        # gamma_val = line_inv[:, :Upgrade_Lvl]
        # storage_inv = CSV.read(storage_file, DataFrame)

        y_val = [gamma_val, s_energy_int_val]
        export_investments_csv(data, gamma_val, s_energy_int_val, output_dir=joinpath(simdir,"benders_output"), file_suffix="$iter")

        # Solve subproblem with fixed investments
        set_sub_investments!(sub, y_val)
        # Phase 1: check feasibility of no load shed
        set_underserved_objective!(sub) 
        solve_subproblem_ptdf!(sub; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6, logging="Eval Point iter $iter")
        total_ue = objective_value(sub)
        duals = extract_subproblem_duals(sub)

        # If failed Phase 1, infeasible due to load shed
        if total_ue > 1e-2
            println("[DEBUG] Master iter $iter: add FEASIBILITY cut, load shed detected")
            export_duals_csv(simdir, data, duals, iter)

            # Save progress to CSV
            filename = joinpath(simdir, "output", "benders_progress.csv")
            benders_ptdf_write_to_csv(filename, y_val, master_obj, theta_val, total_ue * 365 * 5e6, total_ue)
            add_benders_cut_ptdf(master, ue_sum, duals, y, y_val, total_ue)

        else # Passed Phase 1, move onto Phase 2 to check operational objective
            println("[DEBUG] Master iter $iter: add OPTIMALITY cut, there is no load shed")
            set_operational_objective!(sub)
            solve_subproblem_ptdf!(sub; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6)
            phi_val = objective_value(sub)
            duals = extract_subproblem_duals(sub)
            total_ue = vec(sum(value.(sub[:ue]),dims=[1,2,3]))[1]

            export_duals_csv(simdir, data, duals, iter)
            # Save progress to CSV
            filename = joinpath(simdir, "output", "benders_progress.csv")
            benders_ptdf_write_to_csv(filename, y_val, master_obj, theta_val, phi_val, total_ue)
            add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)

            # Update serious step (move trust region)
            push!(master.ext[:y_trust], y_val)

            # Check convergence
            gap = abs(master_obj + phi_val - theta_val - master_obj) / (1e-10 + abs(master_obj))
            push!(master.ext[:gap], gap)
            println("[DEBUG] Benders iteration $(iter+1): Master objective = $(master_obj), theta = $(theta_val), phi = $(phi_val), gap = $(gap)")
        end

        # Update with serious step of the trust/anchor point if satisfy conditions
        if master.ext[:stabilization] == "boxstep"
            if total_ue <= 1e-2 || total_ue < master.ext[:total_ue][end]
                push!(master.ext[:y_trust], y_val)
                push!(master.ext[:total_ue], total_ue)
            end
        elseif master.ext[:stabilization] == "trust_region"
            if total_ue > 1e-2
                new_l1_radius = master.ext[:l1_radius][end] * 2
                println("[DEBUG] Load shed $total_ue detected, expanding transmission l1_radius to $(new_l1_radius)")
                push!(master.ext[:l1_radius], new_l1_radius)
            else
                current_total_obj = master_obj + phi_val - theta_val
                if (current_total_obj < master.ext[:total_obj][end])
                    println("[DEBUG] Cheaper total objective $(current_total_obj) detected, taking serious step")
                    push!(master.ext[:y_trust], y_val)
                    push!(master.ext[:total_obj], current_total_obj)

                    new_l1_radius = maximum([1, master.ext[:l1_radius][end] / 2])
                    println("[DEBUG] Serious step, so reducing l1_radius to $(new_l1_radius)")
                    push!(master.ext[:l1_radius], new_l1_radius)
                else
                    new_l1_radius = master.ext[:l1_radius][end] * 2
                    println("[DEBUG] No load shed or improvement in objective, expanding transmission l1_radius to $(new_l1_radius)")
                    push!(master.ext[:l1_radius], new_l1_radius)
                end
            end
        end

        # Update iter count and save master problem with metadata/extensions
        master.ext[:iter] += 1
        iter = master.ext[:iter]
        # save_master_problem(master, simdir)      
    end

    if !converged
        println("[DEBUG] Benders decomposition did not converge after $max_iterations iterations")
    end
    
    # Return final solution
    return [gamma_val, s_energy_int_val]
end

function solve_pareto_and_cut(master, sub_model, data, theta, y, y_eval, y_core)
    dual_sub_model = make_dual_subproblem(sub_model, data, y_eval)
                
    # modify dual to instead solve for pareto optimal duals
    optimal_obj = objective_bound(sub_model)
    dual_sub_model = solve_for_pareto(dual_sub_model, sub_model, data, optimal_obj, y_eval, y_core)

    # add cut to the master
    add_pareto_cut(master, dual_sub_model, sub_model, data, theta, y)
end

function benders_ptdf_solve(simdir)
    setup_simdir(simdir)

    data = set_up_data(simdir)
    # Create master problem (or load if it exists)
    master, y, theta = load_master_problem(simdir, data)
    
    # Get optimal solution through Benders iterations
    gamma_val, s_energy_int_val = benders_iteration_ptdf(simdir, master, y, theta, data)
    
    # Save final solution
    save_investment_results(simdir, gamma_val, s_energy_val * data["param"]["storage_energy_size"])
    
    return objective_value(master)
end

function save_investment_results(simdir, gamma_val, s_energy_val)
    # Save line investments
    line_df = DataFrame(Upgrade_Lvl = gamma_val)
    CSV.write(joinpath(simdir, "output", "line_investments.csv"), line_df)
    
    # Save storage investments
    storage_df = DataFrame(
        Storage_Energy = s_energy_val
    )
    CSV.write(joinpath(simdir, "output", "storage_investments.csv"), storage_df)
end

function benders_ptdf_write_to_csv(filename, y_val, master_obj, theta_val, phi_val, total_ue)
    # benders_ptdf_write_to_csv(filename, y_val, y_trust_val, master_obj, theta_val, master_obj_trust, theta_trust_val, phi_val)
    # Decompose y_val
    gamma_val, s_energy_int_val = y_val

    # Prepare the base data row
    data_dict = OrderedDict(
        "master_obj_lb" => master_obj,
        "theta_val" => theta_val,
        "phi_val" => phi_val,
        "master_obj_ub" => master_obj + phi_val - theta_val,
        "total_ue" => total_ue,
        "total_line_upgrades" => count(x -> x >=0.015, gamma_val),
        "sum_line_upgrades" => sum(gamma_val),
        "total_storage_energy" => sum(s_energy_int_val),
        "total_storage_count" => count(x -> x >=0.015, s_energy_int_val)
    )
    
    # Create DataFrame from dictionary
    df = DataFrame([data_dict])
    
    # Check if file exists to determine if we need headers
    file_exists = isfile(filename)
    
    if file_exists
        # Read existing headers to ensure compatibility
        existing_headers = names(CSV.read(filename, DataFrame; limit=1))
        
        # Write to CSV, appending to it
        CSV.write(filename, df; append=true, header=false)
    else
        # New file, write with headers
        CSV.write(filename, df)
    end
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

function export_duals_csv(simdir, data, duals, iter)
    # Create the output directory for duals
    mkpath(joinpath(simdir, "output", "duals"))
    
    # Extract transmission duals and save to CSV
    transmission_duals = duals[1]
    transmission_df = DataFrame(transmission_duals = transmission_duals)
    CSV.write(joinpath(simdir, "output", "duals", "transmission_duals_$iter.csv"), transmission_df)
    
    # Extract storage duals and save to CSV
    storage_duals = duals[2]
    storage_df = DataFrame(storage_duals = storage_duals)
    CSV.write(joinpath(simdir, "output", "duals", "storage_duals_$iter.csv"), storage_df)
end

function get_lmp_discharge_value(sub_model)
    tracked_constraints = sub_model.ext[:tracked_constraints]
    ptdf = sub_model.ext[:PTDF]
    E, N = size(ptdf)

    ptdf_constraints = sub_model[:ptdf_flow]
    ptdf_duals = Dict{String, Float64}()
    for (constraint_name, constraint_ref) in ptdf_constraints
        dual_value = dual(constraint_ref)
        ptdf_duals[constraint_name] = dual_value
    end

    # dual_pb = [dual(sub_model[:power_balance][i]) for i=1:24]

    nodal_lmp = zeros(Float64, N, 24)
    for (arc, rep, hour, ub) in keys(tracked_constraints)
        line_dual = nothing
        if ub == 1
            line_dual = ptdf_duals["$(arc)_$(rep)_$(hour)_ub"] * -1.0
        else
            line_dual = ptdf_duals["$(arc)_$(rep)_$(hour)_lb"]
        end
        
        for node in 1:N
            multiplier = ptdf[arc, node]
            if multiplier < 0.0 || multiplier >= 0.05
                # Matrix indexing: nodal_lmp[node, hour]
                nodal_lmp[node, hour] += ptdf[arc, node] * line_dual
            end
        end
    end

    # Add global power balance dual contribution
    """
    for node in 1:N
        for hour in 1:24
            nodal_lmp[node, hour] += dual_pb[hour]
        end
    end
    """

    nodal_lmp *= 0.25

    nodal_lmp_collapsed = sum(nodal_lmp, dims=2)  # Sum across columns (hours)

    positive_row_indices = [idx[1] for idx in findall(nodal_lmp_collapsed .> 0)]
    positive_values = nodal_lmp_collapsed[positive_row_indices]
    large_value_indices = [idx[1] for idx in findall(nodal_lmp_collapsed .> 1e8)]

    result = zeros(N)
    result[positive_row_indices] = -1.0 * positive_values
    return result
end

function create_arcs_to(data)
    arcs_to = Dict{String, Vector{Int}}()
    
    # Iterate through all branches/arcs
    for (arc_idx_str, branch_data) in data["branch"]
        arc_idx = parse(Int, arc_idx_str)
        from_bus = branch_data["f_bus"]
        from_bus_str = string(from_bus)
        
        # Add this arc to the from_bus's list
        if haskey(arcs_to, from_bus_str)
            push!(arcs_to[from_bus_str], arc_idx)
        else
            arcs_to[from_bus_str] = [arc_idx]
        end
    end
    
    return arcs_to
end

function add_storage_minimum!(master)
    data = master.ext[:data]
    storage_file = joinpath(data["param"]["core_point_simdir"], "output", "storage_investments.csv")
    storage_inv = CSV.read(storage_file, DataFrame)
    s_energy_cont_val = storage_inv[:, :Storage_Energy]    
    total_s_energy_cont = sum(s_energy_cont_val)
    total_s_energy_int = ceil.(total_s_energy_cont / data["param"]["storage_energy_size"])

    s_energy_int = master[:s_energy_int]
    N = length(data["bus"])

    @constraint(master,
        total_storage,
        sum(s_energy_int[i] for i in 1:N) >= total_s_energy_int)
end

function add_trust_region!(master)
    if master.ext[:stabilization] != "trust_region"
        return
    end

    data = master.ext[:data]
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    gamma, s_energy_int = master[:gamma], master[:s_energy_int]
    abs_diff, trans_abs_diff = master[:abs_diff], master[:trans_abs_diff]
    
    if master.ext[:iter] == 0
        # warm start initially to the exact over invested point
        gamma_core, s_energy_int_core = zeros(E), zeros(N)
        if master.ext[:warmstart]
            gamma_core, s_energy_int_core = get_over_invested_point(data["param"]["core_point_simdir"], data)
        end
        y_core = [gamma_core, s_energy_int_core]
        push!(master.ext[:y_trust], y_core)
        master.ext[:l1_radius] = [1]
        
        y_trust = master.ext[:y_trust][end]
        gamma_trust, s_energy_int_trust = y_trust

        @constraint(master,
            trans_abs_diff_ub[a in 1:E],
            trans_abs_diff[a] >= gamma[a] - gamma_trust[a])
        @constraint(master,
            trans_abs_diff_lb[a in 1:E],
            trans_abs_diff[a] >= gamma_trust[a] - gamma[a])

        @constraint(master,
            abs_diff_ub[i in 1:N],
            abs_diff[i] >= s_energy_int[i] - s_energy_int_trust[i])
        @constraint(master,
            abs_diff_lb[i in 1:N],
            abs_diff[i] >= s_energy_int_trust[i] - s_energy_int[i])

        @constraint(master,
            trans_abs_diff_total,
            sum(trans_abs_diff[a] for a in 1:E) == 0)
        @constraint(master,
            abs_diff_total,
            sum(abs_diff[i] for i in 1:N) == 0)
        return
    else
        remove_trust_region!(master)
    end
    
    # Move these assignments outside the if-else block to ensure they're always defined
    y_trust = master.ext[:y_trust][end]
    gamma_trust, s_energy_int_trust = y_trust

    @constraint(master,
        trans_abs_diff_ub[a in 1:E],
        trans_abs_diff[a] >= gamma[a] - gamma_trust[a])
    @constraint(master,
        trans_abs_diff_lb[a in 1:E],
        trans_abs_diff[a] >= gamma_trust[a] - gamma[a])
    
    @constraint(master,
        abs_diff_ub[i in 1:N],
        abs_diff[i] >= s_energy_int[i] - s_energy_int_trust[i])
    @constraint(master,
        abs_diff_lb[i in 1:N],
        abs_diff[i] >= s_energy_int_trust[i] - s_energy_int[i])

    @constraint(master,
            trans_abs_diff_total,
            sum(trans_abs_diff[a] for a in 1:E) <= master.ext[:l1_radius][end])
    @constraint(master,
            abs_diff_total,
            sum(abs_diff[i] for i in 1:N) == 2)
end

function add_box_step!(master, box_size)
    if master.ext[:stabilization] != "boxstep"
        return
    end

    data = master.ext[:data]
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    if master.ext[:iter] == 0
        box_size = 0
        gamma_core, s_energy_int_core = zeros(E), zeros(N)
        
        if master.ext[:warmstart]
            gamma_core, s_energy_int_core = get_over_invested_point(data["param"]["core_point_simdir"], data)
        end
        y_core = [gamma_core, s_energy_int_core]
        push!(master.ext[:y_trust], y_core)
    else
        remove_trust_region!(master)
    end

    y_trust = master.ext[:y_trust][end]
    gamma_trust, s_energy_int_trust = y_trust
    gamma, s_energy_int = master[:gamma], master[:s_energy_int]

    @constraint(master,
        trans_box_ub[a in 1:E],
        gamma[a] - gamma_trust[a] <= box_size)

    @constraint(master,
        trans_box_lb[a in 1:E],
        gamma_trust[a] - gamma[a] <= box_size)

    @constraint(master,
        stor_box_ub[i in 1:N],
        s_energy_int[i] - s_energy_int_trust[i] <= box_size)

    @constraint(master,
        stor_box_lb[i in 1:N],
        s_energy_int_trust[i] - s_energy_int[i] <= box_size)
end

function remove_trust_region!(master)
    stab_method = master.ext[:stabilization]

    if stab_method == "trust_region"
        constraint_names = [
            :trans_abs_diff_ub, :trans_abs_diff_lb, :abs_diff_ub, :abs_diff_lb, :abs_diff_total, :trans_abs_diff_total
        ]
    elseif stab_method == "boxstep"
        constraint_names = [
        :trans_box_ub, :trans_box_lb, :stor_box_lb, :stor_box_ub
        ]
    else
        return
    end

    obj_dict = object_dictionary(master)

    for name in constraint_names
        if haskey(obj_dict, name)
            try
                constraint_obj = obj_dict[name]

                # JuMP.delete works on both single constraints and arrays
                JuMP.delete(master, constraint_obj)
                JuMP.unregister(master, name)
                println("Removed constraint: $name")

            catch e
                println("Warning: Could not remove $name: $e")
                # Force unregister
                try
                    JuMP.unregister(master, name)
                catch
                end
            end
        else
            println("Constraint $name not found in model")
        end
    end
end

function get_max_flow_capacity(sub, data)
    E = length(data["branch"])
    K = data["param"]["num_cap_upgrades_max"]
    gamma = value.(sub[:gamma])

    max_flow_capacity = zeros(E)

    for a in 1:E
        cap_increment = get_capacity_increment(data, a)
        line_limit_base = data["branch"]["$a"]["rate_a"]
        max_flow_capacity[a] = line_limit_base + K * cap_increment
        # max_flow_capacity[a] = line_limit_base + gamma[a] * cap_increment
    end

    return max_flow_capacity
end

function compute_flow_overhead(sub, data, hour)
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])
    nodal_injections = zeros(N)
    gen_bus_map = sub.ext[:gen_bus_map]
    pg = value.(sub[:pg])
    ue = value.(sub[:ue])
    ch, dis = value.(sub[:ch]), value.(sub[:dis])

    for i in 1:N
        nodal_injections[i] = (sum(pg[1, g, hour] for g in 1:G if gen_bus_map[g] == i; init=0.0)
                - data["bus"]["$i"]["load"]["1"][hour]
                + ue[1, i, hour] - ch[1, i, hour] + dis[1, i, hour])
    end

    ptdf = sub.ext[:PTDF]
    flows = ptdf * nodal_injections

    max_flow_capacity = get_max_flow_capacity(sub, data)

    max_inj_before_violation = zeros(N)

    m_minus_f = max_flow_capacity - flows
    minus_m_minus_f = - max_flow_capacity - flows

    for i in 1:N
        injections_allowed = zeros(E)
        for a in 1:E
            injections_allowed[a] = ptdf[a,i] > 0 ? m_minus_f[a] / ptdf[a,i] : (ptdf[a,i] < 0 ? minus_m_minus_f[a] / ptdf[a,i] : Inf)
        end

        max_inj_before_violation[i] = minimum(injections_allowed)
    end
end
