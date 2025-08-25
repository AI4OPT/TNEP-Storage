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

    # siting abs diff
    JuMP.@variable(master, abs_diff[1:N] >= 0)

    # sizing abs diff
    JuMP.@variable(master, sizing_abs_diff >= 0)

    # transmission abs diff
    JuMP.@variable(master, trans_abs_diff[a=1:E] >= 0)

    #
    #   II. Constraints
    #    

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

    # energy rating only if storage installed
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy_int[i] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
    )

    #
    #   III. Objective
    #
    # Start with the base objective terms
    obj_expr = sum(s_energy_int[i] * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            theta
    y = gamma, s_energy_int

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)

    return master, y, theta
end

function remove_trust_region!(master)
    constraint_names = [
        :trans_abs_diff_ub, :trans_abs_diff_lb, :transmission_trust_region,
        :abs_diff_ub, :abs_diff_lb, :siting_trust_region,
        :sizing_abs_diff_ub, :sizing_abs_diff_lb, :sizing_trust_region
    ]
    
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

function add_trust_region!(master, data)
    # If first iteration, get the starting point
    if master.ext[:iter] == 0
        # TODO: write new over-invested core function
        gamma_core, s_energy_int_core = get_over_invested_point(data["param"]["core_point_simdir"], data)
        y_core = [gamma_core, s_energy_int_core]
        push!(master.ext[:y_trust], y_core)
        push!(master.ext[:trust_siting], 0.0)
        push!(master.ext[:trust_sizing], 0.0)
        push!(master.ext[:trust_transmission], 0.0)
    end
    
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    
    y_trust = master.ext[:y_trust][end]
    gamma_trust, s_energy_int_trust = y_trust
    gamma, s_energy_int = master[:gamma], master[:s_energy_int]
    abs_diff = master[:abs_diff]
    sizing_abs_diff = master[:sizing_abs_diff]
    trans_abs_diff = master[:trans_abs_diff]
    
    # Transmission trust region
    @constraint(master,
        trans_abs_diff_ub[a in 1:E],
        trans_abs_diff[a] >= gamma[a] - gamma_trust[a])
    
    @constraint(master,
        trans_abs_diff_lb[a in 1:E],
        trans_abs_diff[a] >= gamma_trust[a] - gamma[a])
    
    @constraint(master, 
        transmission_trust_region, 
        sum(trans_abs_diff[a] for a in 1:E) <= master.ext[:trust_transmission][end])
    
    # Sizing trust region (Note: your comment says "sizing" but constraints use "siting" variables)
    @constraint(master,
        abs_diff_ub[i in 1:N],
        abs_diff[i] >= s_energy_int[i] - s_energy_int_trust[i])
    
    @constraint(master,
        abs_diff_lb[i in 1:N],
        abs_diff[i] >= s_energy_int_trust[i] - s_energy_int[i])
    
    @constraint(master, 
        siting_trust_region, 
        sum(abs_diff[i] for i in 1:N) <= master.ext[:trust_siting][end])
    
    # Siting trust region (Note: your comment says "siting" but uses sizing variables)
    @constraint(master,
        sizing_abs_diff_ub,
        sizing_abs_diff >= sum(s_energy_int[i] for i in 1:N) - sum(s_energy_int_trust[i] for i in 1:N))
    
    @constraint(master,
        sizing_abs_diff_lb,
        sizing_abs_diff >= sum(s_energy_int_trust[i] for i in 1:N) - sum(s_energy_int[i] for i in 1:N))
    
    @constraint(master, 
        sizing_trust_region, 
        sizing_abs_diff <= master.ext[:trust_sizing][end])
    
    # Verify all constraints were added successfully
    println("Trust region constraints added successfully for iteration ", master.ext[:iter])
end

function solve_subproblem_ptdf(simdir, y_val, data, tracked_constraints; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6, logging=nothing)

    if !isnothing(logging)
        println("[DEBUG] $logging: Beginning subproblem solve")
    end

    # unpack investment decisions
    gamma_val, s_energy_int_val = y_val

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
    @variable(sub, s_energy_int[i=1:N])

    #
    #   II. Constraints
    #

    # FIX INVESTMENT DECISIONS FROM MASTER PROBLEM (DUALS USED FOR CUTS)
    @constraint(sub, master_gamma[a=1:E], gamma[a] == gamma_val[a])
    @constraint(sub, master_energy[i=1:N], s_energy_int[i] == s_energy_int_val[i])

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
                flow_val = flows[a, r, t]
                violation_amount = abs(flow_val) - line_limit_fixed
                ub = (flow_val >= 0)

                if get(sub.ext[:tracked_constraints], (a, r, t, ub), false)
                    continue
                end
        
                if violation_amount > sub.ext[:solve_metadata][:ptdf_tol]
                    push!(violations, (a, r, t, ub, violation_amount))
                end
            end
        end

        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        n_added = 0

        # Add up to the limit
        max_to_add = sub.ext[:solve_metadata][:max_ptdf_per_iteration]
        for (a, r, t, ub, _) in Iterators.take(sorted_violations, max_to_add)

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
        dual_energy = [dual(master_energy[i]) for i=1:N]
        
        duals = (dual_gamma, dual_energy)

        return sub, objective_value(sub), duals, sub.ext[:tracked_constraints]
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

    while !converged && iter < max_iterations
        # Solve master problem
        println("[DEBUG] Solving master problem: iteration $iter")
        if iter > 0
            remove_trust_region!(master)
        end
        optimize!(master)
        gamma_no_trust_val = value.(master[:gamma])
        s_energy_int_no_trust_val = value.(master[:s_energy_int])
        y_no_trust = [gamma_no_trust_val, s_energy_int_no_trust_val]
        export_investments_csv(data, gamma_no_trust_val, s_energy_int_no_trust_val, output_dir=joinpath(simdir,"benders_output"), file_suffix="$iter")

        # Save lower bound (no trust region)
        master_obj_lb = objective_value(master) 

        # Add trust region
        add_trust_region!(master, data)
        optimize!(master)
        master_obj_trust = objective_value(master)
        
        # Get master solution
        gamma_val = value.(master[:gamma])
        s_energy_int_val = value.(master[:s_energy_int])
        theta_val = value(theta)
        y_raw = [gamma_val, s_energy_int_val]
        export_investments_csv(data, gamma_val, s_energy_int_val, output_dir=joinpath(simdir,"benders_output"), file_suffix="trust_$iter")

        # Solve subproblem with fixed investments
        sub_model, phi_val, duals, raw_tracked_constraints = solve_subproblem_ptdf(simdir, y_raw, data, tracked_constraints, logging="Eval Point iter $iter")
        # TODO: export_duals_csv(duals, iter)

        # Update tracked constraints
        tracked_constraints = raw_tracked_constraints

        # Save progress to CSV
        filename = joinpath(simdir, "output", "benders_progress.csv")
        benders_ptdf_write_to_csv(filename, y_no_trust, y_raw, master_obj_lb, theta_val, phi_val, master_obj_trust + phi_val - theta_val)

        # Check convergence
        gap = abs(master_obj_trust + phi_val - theta_val - master_obj_lb) / (1e-10 + abs(master_obj_lb))
        push!(master.ext[:gap], gap)
        println("[DEBUG] Benders iteration $(iter+1): Master objective = $(objective_value(master)), theta = $(theta_val), phi = $(phi_val), gap = $(gap)")

        # if gap < tolerance
        if false
            return
        else
            # Add new Benders cut to master problem
            if raw_model !== nothing
                add_appropriate_cut(master, raw_model, data, theta, y, y_raw, y_core, raw_duals, raw_phi_val)
            end
            add_appropriate_cut(master, sub_model, data, theta, y, y_raw, y_raw, duals, phi_val)
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

function benders_ptdf_write_to_csv(filename, y_og, y_val, master_obj_lb, theta_val, phi_val, master_obj_trust)
    # benders_ptdf_write_to_csv(filename, master_obj_lb, theta_val, phi_val, master_obj_trust)
    # Decompose y_val
    gamma_og, s_energy_int_og = y_og
    gamma_val, s_energy_int_val = y_val

    # Prepare the base data row
    data_dict = OrderedDict(
        "master_obj_lb" => master_obj_lb,
        "theta_val" => theta_val,
        "phi_val" => phi_val,
        "master_obj_trust" => master_obj_trust,
        "total_line_upgrades" => count(x -> x >=0.015, gamma_og),
        "sum_line_upgrades" => sum(gamma_og),
        "total_storage_energy" => sum(s_energy_int_og),
        "total_storage_count" => count(x -> x >=0.015, s_energy_int_og),
        "total_line_upgrades_trust" => count(x -> x >=0.015, gamma_val),
        "sum_line_upgrades_trust" => sum(gamma_val),
        "total_storage_energy_trust" => sum(s_energy_int_val),
        "total_storage_count_trust" => count(x -> x >=0.015, s_energy_int_val)
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