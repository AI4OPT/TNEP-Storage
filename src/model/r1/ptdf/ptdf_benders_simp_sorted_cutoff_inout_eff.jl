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
include("ptdf_benders_transport_flow_subproblem.jl")
include("ptdf_core_point.jl")
include("ptdf_benders_pareto.jl")

function define_master_ptdf(data::Dict{String, Any})
    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)
    set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Create extension for tracking iterations of y_raw, y_eval, and y_core
    master.ext[:y_raw] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_eval] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_core] = Vector{Vector{Vector{Float64}}}()

    # Create extension for tracking whether transport flow should be embedded
    master.ext[:transport] = haskey(data["param"], "embed_transport_flow") && data["param"]["embed_transport_flow"] == true

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
    if haskey(data["param"], "relaxed_first_stage") && data["param"]["relaxed_first_stage"] == true
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
    JuMP.@variable(master, theta >= 0)

    #
    #   II. Constraints
    #    

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
    # Is set in the Bender's loop
    y = gamma, s_power, s_energy

    return master, y, theta

end

function update_master_objective!(master, data, y, theta, iter)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    #
    #   III. Objective
    #
    gamma, s_power, s_energy = y

    # Start with the base objective terms
    obj_expr = sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            theta

    # L2 (quadratic) regularization
    if haskey(data["param"], "reg_penalty")
        gamma_reg, s_power_reg, s_energy_reg = compute_regularization_point(master, data, iter)
        obj_expr += data["param"]["reg_penalty"] * sum((s_power[i] - s_power_reg[i])^2 for i in 1:N)
        if haskey(data["param"], "trans_reg_penalty")
            obj_expr += data["param"]["trans_reg_penalty"] * sum((gamma[a] - gamma_reg[a])^2 for a in 1:E)
        end
    end

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)
end

function solve_subproblem_ptdf(simdir, y_val, data, tracked_constraints; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6, logging=nothing)

    if !isnothing(logging)
        println("[DEBUG] $logging: Beginning subproblem solve")
    end

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
        soc[r,i,1] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] + ch[r,i,1] * data["param"]["bess_efficiency"] - dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    @constraint(sub,
        soc_end[r=1:R, i=1:N],
        soc[r,i,T] == get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i]
    )

    # SOC limits
    @constraint(sub, 
        soc_energy_ub[r=1:R, i=1:N, t=1:T],
        soc[r,i,t] <= s_energy[i]
    )
    
    # Charging/discharging limits
    @constraint(sub,
        charge_discharge_lb[r=1:R, i=1:N, t=1:T],
        dis[r,i,t] <= s_power[i]
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
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i=1:N) + 
                    sum(get(data["param"], "storage_operation_cost", 0.0) * (ch[r,i,t] + dis[r,i,t]) for i=1:N)
                for t=1:T)
            )
        for r=1:R) * operational_weight
    )

    # Apply warm-start using tracked constraints
    if !isempty(tracked_constraints)
        println("[DEBUG] Beginning optimized warm-start for subproblem")
        
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

            # Pre-compute the line limits for all constraints in this group
            line_limit = line_limit_base + sub[:gamma][a] * cap_increment
            
            # Generate constraints efficiently
            for (r, t) in rep_time_pairs
                # Create flow expression once and reuse
                flow_expr = sum(ptdf_row[i] * (
                    sum(sub[:pg][r,g,t] for g in 1:G if gen_bus_map[g] == i; init=0.0) -
                    data["bus"]["$i"]["load"]["$r"][t] +
                    sub[:ue][r,i,t] -
                    sub[:ch][r,i,t] +
                    sub[:dis][r,i,t]
                ) for i in 1:N)
                
                # Add both constraints
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit)
                sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit)
            end
        end
        
        println("[DEBUG] Warm-started subproblem with $total_constraints constraints")
    end

    # Solve using lazy PTDF constraints
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < sub.ext[:solve_metadata][:max_ptdf_iterations]
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
                if get(sub.ext[:tracked_constraints], (a, r, t), false)
                    continue
                end
        
                flow_val = flows[a, r, t]
                violation_amount = abs(flow_val) - line_limit_fixed
        
                if violation_amount > sub.ext[:solve_metadata][:ptdf_tol]
                    push!(violations, (a, r, t, violation_amount))
                end
            end
        end

        sorted_violations = sort(violations, by = x -> -x[4])
        n_violated = length(sorted_violations)
        n_added = 0

        # Add up to the limit
        max_to_add = sub.ext[:solve_metadata][:max_ptdf_per_iteration]
        for (a, r, t, _) in Iterators.take(sorted_violations, max_to_add)

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
            sub[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(sub, flow_expr <= line_limit_expr)
            sub[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(sub, flow_expr >= -line_limit_expr)
        
            sub.ext[:tracked_constraints][(a, r, t)] = true
            n_added += 1
        end

        if n_added > 0
            save_tracked_constraints(simdir, sub, n_violated)
        end
        solved = (n_violated == 0)
        niter += 1

        # println("approx. progress $(checked_count / E)")
        println("[DEBUG] PTDF Subproblem Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    solve_time = time() - t0
    sub.ext[:solve_metadata][:solve_time] = solve_time

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
        println("[DEBUG] Loading $(nrow(tracked_df)) previously tracked constraints")
        
        for row in eachrow(tracked_df)
            if row.tracked
                tracked_constraints[(row.arc, row.rep, row.time)] = true
            end
        end
    end

    while !converged && iter < max_iterations
        # Solve master problem
        println("[DEBUG] Solving master problem: iteration $iter")
        update_master_objective!(master, data, y, theta, iter)
        optimize!(master)
        
        if termination_status(master) != MOI.OPTIMAL
            error("Master problem failed to solve optimally: $(termination_status(master))")
        end
        
        # Get master solution
        gamma_val = value.(gamma)
        s_power_val = value.(s_power)
        s_energy_val = value.(s_energy)
        theta_val = value(theta)
        y_raw = [gamma_val, s_power_val, s_energy_val]
        push!(master.ext[:y_raw], y_raw)
        export_investments_csv(data, gamma_val, s_power_val, s_energy_val, output_dir=joinpath(simdir,"benders_output"), file_suffix="$iter")

        y_eval, y_core = compute_eval_core_points(master, data, iter)


        # Solve core point subproblem to find the PTDF constraints
        core_model, core_phi_val, core_duals, core_tracked_constraints = solve_subproblem_ptdf(simdir, y_core, data, tracked_constraints, logging="Core Point iter $iter")
        # Solve subproblem with fixed investments
        sub_model, phi_val, duals, eval_tracked_constraints = solve_subproblem_ptdf(simdir, y_eval, data, core_tracked_constraints, logging="Eval Point iter $iter")

        # Update tracked constraints
        tracked_constraints = eval_tracked_constraints

        # Save progress to CSV
        filename = joinpath(simdir, "output", "benders_progress.csv")
        # The benders log will show the perturbed / stabilized point
        benders_ptdf_write_to_csv(filename, objective_value(master), theta_val, phi_val, y_eval)

        # Check convergence
        gap = abs(theta_val - phi_val) / (1e-10 + abs(phi_val))
        println("[DEBUG] Benders iteration $(iter+1): Master objective = $(objective_value(master)), theta = $(theta_val), phi = $(phi_val), gap = $(gap)")

        # if gap < tolerance
        if false
            # println("[DEBUG] Convergence detected with stabilized solution. Performing final cut check at discrete solution...")
            println("[DEBUG] Convergence detected...")

            """
            # Add new Benders cut to master problem
            add_appropriate_cut(master, sub_model, data, theta, y, y_eval, y_core, duals, phi_val)
            optimize!(master)
        
            # Get updated unstabilized master solution
            gamma_val = value.(gamma)
            s_power_val = value.(s_power)
            s_energy_val = value.(s_energy)
            theta_val = value(theta)
        
            # Final subproblem solve at discrete solution
            y_final = [gamma_val, s_power_val, s_energy_val]
            _, final_phi_val, final_duals, final_tracked_constraints = solve_subproblem_ptdf(simdir, y_final, data, tracked_constraints)
            tracked_constraints = final_tracked_constraints
        
            # Save progress to CSV
            filename = joinpath(simdir, "output", "benders_progress.csv")
            benders_ptdf_write_to_csv(filename, objective_value(master), theta_val, final_phi_val, y_final)
        
            final_gap = abs(theta_val - final_phi_val) / (1e-10 + abs(final_phi_val))"""
            final_gap = gap
        
            if final_gap < tolerance
                println("[DEBUG] Final gap at discrete solution = $(final_gap) < $(tolerance). Converged after $(iter+1) iterations.")
                converged = true
            else
                # If still not converged, add another cut and continue
                add_appropriate_cut(master, sub_model, data, theta, y, y_eval, y_core, duals, phi_val)
            end
        else
            # Add new Benders cut to master problem
            add_appropriate_cut(master, sub_model, data, theta, y, y_eval, y_core, duals, phi_val)
        end
        
        iter += 1
    end

    if !converged
        println("[DEBUG] Benders decomposition did not converge after $max_iterations iterations")
    end
    
    # Return final solution
    return [gamma_val, s_power_val, s_energy_val]
end

# Helper function to add the appropriate cut
function add_appropriate_cut(master, sub_model, data, theta, y, y_eval, y_core, duals, phi_val)
    if haskey(data["param"], "benders_technique") && data["param"]["benders_technique"] == "pareto"
        solve_pareto_and_cut(master, sub_model, data, theta, y, y_eval, y_core)
    elseif haskey(data["param"], "benders_technique") && data["param"]["benders_technique"] == "both"
        solve_pareto_and_cut(master, sub_model, data, theta, y, y_eval, y_core)
        add_benders_cut_ptdf(master, theta, duals, y, y_eval, phi_val)
    else
        add_benders_cut_ptdf(master, theta, duals, y, y_eval, phi_val)
    end
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

function benders_ptdf_write_to_csv(filename, master_obj, theta_val, phi_val, y_val, phirel_val=nothing)
    # Decompose y_val
    gamma_val, s_power_val, s_energy_val = y_val

    # Prepare the base data row
    data_dict = OrderedDict(
        "master_objective" => master_obj,
        "theta_val" => theta_val,
        "phi_val" => phi_val,
        "total_line_upgrades" => count(x -> x > 0, gamma_val),
        "total_storage_power" => sum(s_power_val),
        "total_storage_energy" => sum(s_energy_val),
        "total_storage_count" => count(x -> x != 0, s_power_val)
    )
    
    # Add phirel_val if provided
    if phirel_val !== nothing
        data_dict["phirel_val"] = phirel_val
    end
    
    # Create DataFrame from dictionary
    df = DataFrame([data_dict])
    
    # Check if file exists to determine if we need headers
    file_exists = isfile(filename)
    
    if file_exists
        # Read existing headers to ensure compatibility
        existing_headers = names(CSV.read(filename, DataFrame; limit=1))
        
        # Check if phirel_val column exists in the file but not in our data
        if "phirel_val" in existing_headers && !("phirel_val" in keys(data_dict))
            # Add empty phirel_val to maintain column structure
            df.phirel_val = [missing]
        end
        
        # Write to CSV, appending to it
        CSV.write(filename, df; append=true, header=false)
    else
        # New file, write with headers
        CSV.write(filename, df)
    end
end