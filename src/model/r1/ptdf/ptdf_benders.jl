using JuMP
using Gurobi
using CSV, DataFrames
include("../../../helpers/compute_gen_cost.jl")
include("../rate_a_zero.jl")
include("base_ptdf.jl")

function define_master_ptdf(data::Dict{String, Any})
    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    fbus::Vector{Int} = [data["branch"]["$a"]["f_bus"] for a in 1:E]                # from bus
    tbus::Vector{Int} = [data["branch"]["$a"]["t_bus"] for a in 1:E]                # to bus

    #
    #   I. Variables
    #

    # investment level of capacity upgrade
    JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)

    # binary variable for installation of storage
    # JuMP.@variable(master, sigma[i=1:N], Bin)

    # power rating of storage
    JuMP.@variable(master, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(master, s_energy[i=1:N] >= 0)

    # subproblem objective(s)
    JuMP.@variable(master, theta >= 0)

    #
    #   II. Constraints
    #
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(master, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # energy rating only if storage installed
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= data["param"]["max_energy_rating"]
    )

    # power rating only if storage installed
    JuMP.@constraint(master, 
        installed_power_ub[i in 1:N],
        s_power[i] <= data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    JuMP.@constraint(master, 
        short_duration[i in 1:N],
        s_energy[i] == 4.0 * s_power[i]
    )

    #
    #   III. Objective
    #
    JuMP.@objective(master, Min,
    sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) + 
    sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
    theta
    )

    y = gamma, s_power, s_energy

    return master, y, theta

end

function solve_subproblem_ptdf(simdir, y_val, data, tracked_constraints, max_ptdf_iterations=64, max_ptdf_per_iteration=32, ptdf_tol=1e-6)
    # unpack investment decisions
    gamma_val, s_power_val, s_energy_val = y_val

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

    # Set up the PTDF metadata
    sub.ext[:solve_metadata] = Dict(
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
        :ptdf_tol => ptdf_tol
    )
    
    # Compute PTDF matrix
    sub.ext[:PTDF] = do_all_ptdf(data)

    # Initialize tracked constraints with the passed dictionary
    sub.ext[:tracked_constraints] = deepcopy(tracked_constraints)
    
    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    sub.ext[:rate_a_nonzero] = Set(parse(Int, x) for x in rate_a_nonzero)

    # Container for flow constraints
    sub[:ptdf_flow] = Dict{String, ConstraintRef}()

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

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
    @variable(sub, ch[r=1:R, i=1:N, t=1:T])  # Charging/discharging
    @variable(sub, soc[r=1:R, i=1:N, t=1:T] >= 0)  # State of charge

    # Fix investment variables based on master problem decisions
    @variable(sub, gamma[a=1:E])
    @variable(sub, s_power[i=1:N])
    @variable(sub, s_energy[i=1:N])

    #
    #   II. Constraints
    #

    # FIX INVESTMENT DECISIONS FROM MASTER PROBLEM (DUALS USED FOR CUTS)
    @constraint(sub, master_gamma[a=1:E], gamma[a] == gamma_val[a])
    @constraint(sub, master_power[i=1:N], s_power[i] == s_power_val[i])
    @constraint(sub, master_energy[i=1:N], s_energy[i] == s_energy_val[i])

    # Global power balance
    @constraint(sub, 
        power_balance[r=1:R, t=1:T],
        sum(pg[r,g,t] for g=1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i=1:N)
        + sum(ue[r,i,t] for i=1:N)
        - sum(ch[r,i,t] for i=1:N)
        == 0
    )

    # Storage constraints

    # SOC evolution
    @constraint(sub, 
        soc_over_time[r=1:R, i=1:N, t=2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t]
    )
    
    # Initial and final SOC
    @constraint(sub,
        soc_start[r=1:R, i=1:N],
        soc[r,i,1] == 0.5 * s_energy[i] + ch[r,i,1]
    )
    @constraint(sub,
        soc_end[r=1:R, i=1:N],
        soc[r,i,T] == 0.5 * s_energy[i]
    )

    # SOC limits
    @constraint(sub, 
        soc_energy_ub[r=1:R, i=1:N, t=1:T],
        soc[r,i,t] <= s_energy[i]
    )
    
    # Charging/discharging limits
    @constraint(sub,
        charge_discharge_lb[r=1:R, i=1:N, t=1:T],
        -s_power[i] <= ch[r,i,t]
    )
    @constraint(sub,
        charge_discharge_ub[r=1:R, i=1:N, t=1:T],
        ch[r,i,t] <= s_power[i]
    )

    # Objective function
    operational_weight = get(data["param"], "operational_weight", 1)

    @objective(sub, Min,
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g=1:G) +
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i=1:N)
                for t=1:T)
            )
        for r=1:R) * operational_weight
    )

    # Apply warm-start using tracked constraints
    if !isempty(tracked_constraints)
        println("Beginning optimized warm-start for subproblem")
        
        # Convert dictionary to DataFrame for grouping
        arc_list = [k[1] for k in keys(tracked_constraints)]
        rep_list = [k[2] for k in keys(tracked_constraints)]
        time_list = [k[3] for k in keys(tracked_constraints)]
        tracked_df = DataFrame(arc=arc_list, rep=rep_list, time=time_list)
        
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
            rep_time_pairs = [(row.rep, row.time) for row in constraints_batch]
            
            # Generate constraints efficiently
            for (r, t) in rep_time_pairs
                # Expression using the fixed gamma from master problem
                line_limit = line_limit_base + sub[:gamma][a] * cap_increment
                
                # Create flow expression once and reuse
                flow_expr = sum(ptdf_row[i] * (
                    sum(sub[:pg][r,g,t] for g in 1:G if gen_bus_map[g] == i; init=0.0) -
                    data["bus"]["$i"]["load"]["$r"][t] +
                    sub[:ue][r,i,t] -
                    sub[:ch][r,i,t]
                ) for i in 1:N)
                
                # Add both constraints
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit)
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit)
            end
        end
        
        println("Warm-started subproblem with $total_constraints constraints")
    end

    # Solve using lazy PTDF constraints
    solved = false
    niter = 0

    while !solved && niter < sub.ext[:solve_metadata][:max_ptdf_iterations]
        # Solve model
        optimize!(sub)
        
        # Exit if not solved optimally
        st = termination_status(sub)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(sub[:pg])
        ue_values = value.(sub[:ue])
        ch_values = value.(sub[:ch])
        gamma_values = value.(sub[:gamma])

        # Compute all flows at once
        flows = zeros(length(sub.ext[:rate_a_nonzero]), R, T)
        compute_flows!(flows, pg_values, ue_values, ch_values, data, sub.ext[:PTDF])

        # Hybrid approach: Branch-by-branch with vectorized (r,t) checks
        n_violated = 0
        n_added = 0
        checked_count = 0

        for a in sub.ext[:rate_a_nonzero]
            checked_count += 1

            # Pre-compute for this branch
            ptdf_row = sub.ext[:PTDF][a,:]
            cap_increment = get_capacity_increment(data, a)
            line_limit_base = data["branch"]["$a"]["rate_a"]
            line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment
            line_limit_expr = line_limit_base + sub[:gamma][a] * cap_increment

            # Create a matrix of all (r,t) combinations for this branch
            rt_pairs = [(r, t) for r in 1:R for t in 1:T]
            
            # Filter out already tracked combinations (vectorized)
            untracked_rt = filter(rt -> !get(sub.ext[:tracked_constraints], (a, rt[1], rt[2]), false), rt_pairs)

            # Check all flows for this branch at once (vectorized)
            branch_violations = filter(rt -> begin
                r, t = rt
                return abs(flows[a,r,t]) > line_limit_fixed + sub.ext[:solve_metadata][:ptdf_tol]
            end, untracked_rt)

            # Update total violations count
            branch_violation_count = length(branch_violations)
            n_violated += branch_violation_count

            # If any violations, add constraints up to the limit
            if branch_violation_count > 0
                # Determine how many violations to process for this branch
                remaining_slots = sub.ext[:solve_metadata][:max_ptdf_per_iteration] - n_added
                to_process = branch_violations[1:min(length(branch_violations), remaining_slots)]
                
                # Process violations for this branch
                for (r, t) in to_process
                    # Create flow expression only once
                    flow_expr = sum(ptdf_row[i] * (
                        sum(sub[:pg][r,g,t] for g in 1:G if gen_bus_map[g] == i; init=0.0) 
                        - data["bus"]["$i"]["load"]["$r"][t] 
                        + sub[:ue][r,i,t] - sub[:ch][r,i,t]
                    ) for i in 1:N)
                    
                    # Add both constraints with the same expression
                    sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit_expr)
                    sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit_expr)
                    
                    sub.ext[:tracked_constraints][(a,r,t)] = true
                    n_added += 1
                end
                
                # Check if we've hit the limit
                if n_added >= sub.ext[:solve_metadata][:max_ptdf_per_iteration]
                    break
                end
            end
        end
        
        solved = (n_violated == 0)
        niter += 1
        println("approx. progress $(checked_count / E)")
        println("PTDF Subproblem Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Extract duals for Benders cuts
    # Note: Ensure the model is solved to optimality before extracting duals
    if termination_status(sub) ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        dual_gamma = [dual(master_gamma[a]) for a=1:E]
        dual_power = [dual(master_power[i]) for i=1:N]
        
        duals = (dual_gamma, dual_power)
        return sub, objective_value(sub), duals, sub.ext[:tracked_constraints]
    else
        error("Subproblem failed to solve optimally: $(termination_status(sub))")
    end
end

function add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
    # Unpack y variables and duals
    gamma, s_power, s_energy = y
    dual_gamma, dual_power = duals
    gamma_val, s_power_val, s_energy_Val = y_val

    # Add Benders cut
    @constraint(master, 
        theta >= phi_val + 
        sum(dual_gamma[a] * (gamma[a] - gamma_val[a]) for a=1:length(gamma)) +
        sum(dual_power[i] * (s_power[i] - s_power_val[i]) for i=1:length(s_power))
    )
end

function benders_iteration_ptdf(simdir, master, y, theta, data, max_iterations=1000, tolerance=0.01)
    converged = false
    iter = 0
    
    # Unpack y variables
    gamma, s_power, s_energy = y

    # Initialize tracked constraints dictionary to store across iterations
    tracked_constraints = Dict{Tuple{Int,Int,Int}, Bool}()

    # Load previous constraints if available
    if isfile(joinpath(simdir, "tracked_constraints.csv"))
        tracked_df = CSV.read(joinpath(simdir, "tracked_constraints.csv"), DataFrame)
        println("Loading $(nrow(tracked_df)) previously tracked constraints")
        
        for row in eachrow(tracked_df)
            if row.tracked
                tracked_constraints[(row.arc, row.rep, row.time)] = true
            end
        end
    end

    gamma_val, s_power_val, s_energy_val = nothing, nothing, nothing

    while !converged && iter < max_iterations
        # Solve master problem
        optimize!(master)
        
        if termination_status(master) != MOI.OPTIMAL
            error("Master problem failed to solve optimally: $(termination_status(master))")
        end
        
        # Get master solution
        gamma_val = value.(gamma)
        s_power_val = value.(s_power)
        s_energy_val = value.(s_energy)
        theta_val = value(theta)

        y_val = (gamma_val, s_power_val, s_energy_val)

        # Solve subproblem with fixed investments
        sub_model, phi_val, duals, new_tracked_constraints = solve_subproblem_ptdf(simdir, y_val, data, tracked_constraints)

        # Update tracked constraints
        tracked_constraints = new_tracked_constraints

        # Save progress to CSV
        filename = joinpath(simdir, "output", "benders_progress.csv")
        benders_ptdf_write_to_csv(filename, objective_value(master), theta_val, phi_val, gamma_val, s_power_val, s_energy_val)

        # Check convergence
        gap = abs(theta_val - phi_val) / (1e-10 + abs(phi_val))
        println("Benders iteration $(iter+1): Master objective = $(objective_value(master)), theta = $(theta_val), phi = $(phi_val), gap = $(gap)")

        if gap < tolerance
            converged = true
            println("Benders decomposition converged after $(iter+1) iterations")
        else
            # Add new Benders cut to master problem
            add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
        end
        
        iter += 1
    end

    if !converged
        println("Benders decomposition did not converge after $max_iterations iterations")
    end

    save_tracked_constraints(simdir, tracked_constraints)
    
    # Return final solution
    return (gamma_val, s_power_val, s_energy_val)
end

function benders_ptdf_solve(simdir)
    data = JSON.parsefile(joinpath(simdir, "data.json"))
    # Create master problem
    master, y, theta = define_master_ptdf(data)
    
    # Get optimal solution through Benders iterations
    gamma_val, s_power_val, s_energy_val = benders_iteration_ptdf(simdir, master, y, theta, data)
    
    # Save final solution
    save_investment_results(simdir, gamma_val, s_power_val, s_energy_val)
    
    return objective_value(master)
end

function save_investment_results(simdir, gamma_val, s_power_val, s_energy_val)
    # Save line investments
    line_df = DataFrame(Upgrade_Lvl = gamma_val)
    CSV.write(joinpath(simdir, "output", "line_investments.csv"), line_df)
    
    # Save storage investments
    storage_df = DataFrame(
        Storage_Power = s_power_val,
        Storage_Energy = s_energy_val
    )
    CSV.write(joinpath(simdir, "output", "storage_investments.csv"), storage_df)
end

function benders_ptdf_write_to_csv(filename, master_obj, theta_val, phi_val, gamma_val, s_power_val, s_energy_val)
    # Prepare the data row
    df = DataFrame(
        master_objective = [master_obj],
        theta_val = [theta_val],
        phi_val = [phi_val],
        total_line_upgrades = [sum(gamma_val)],
        total_storage_power = [sum(s_power_val)],
        total_storage_energy = [sum(s_energy_val)],
        total_storage_count = [count(x -> x != 0, s_power_val)]
    )

    # Write to CSV, appending to it
    if isfile(filename)
        CSV.write(filename, df; append = true, header = false)
    else
        CSV.write(filename, df)
    end
end

function save_tracked_constraints(simdir, tracked_constraints::Dict{Tuple{Int,Int,Int}, Bool})
    df = DataFrame(
        arc = [k[1] for k in keys(tracked_constraints)],
        rep = [k[2] for k in keys(tracked_constraints)],
        time = [k[3] for k in keys(tracked_constraints)],
        tracked = fill(true, length(tracked_constraints))
    )
    
    CSV.write(joinpath(simdir, "output", "tracked_constraints.csv"), df)
    println("Saved $(nrow(df)) tracked constraints for future warm starts")
end
