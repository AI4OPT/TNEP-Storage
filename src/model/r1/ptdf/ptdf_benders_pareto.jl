function compute_regularization_point(master, data)
    iter = master.ext[:iter]
    params = get(data["param"], "dynamic_lambda", [0.1, 1.0])
    threshold, EXTREME_GAP = params

    # Compute the reg point
    gamma_reg, s_energy_reg = get_rep_day_core_point(data["param"]["core_point_simdir"])

    if !isempty(master.ext[:gap]) # i.e. iter > 0
        @assert length(master.ext[:y_eval]) == length(master.ext[:gap])
        
        # Find the most recent eval point with gap < threshold
        for i in iter:-1:1  # Start from current iteration and go backwards
            if master.ext[:gap][i] < threshold
                y_raw_selected = master.ext[:y_raw][i]
                gamma_reg, s_energy_reg = y_raw_selected
                break
            end
        end
    end

    return gamma_reg, s_energy_reg
end

function get_stabilization_shift(master, data)
    iter = master.ext[:iter]

    # Get lambda parameters
    stabilization_lambda = get(data["param"], "stabilization_lambda", [0.0, 0.0, 0.0])
    starter_lambda = stabilization_lambda[1]
    decrement = stabilization_lambda[2]
    min_lambda = stabilization_lambda[3]
    lambda = nothing
    decr = nothing
    # Get core shift phi parameters
    core_shift_phi = get(data["param"], "core_shift_phi", [0.0, 0.0])
    main_phi = core_shift_phi[1]
    low_phi = core_shift_phi[2]
    phi = nothing

    if isempty(master.ext[:gap]) # i.e. iter == 0
        push!(master.ext[:stabilization_lambda], starter_lambda)
        push!(master.ext[:stabilization_lambda_decr], decrement)
        return starter_lambda, low_phi
    else
        if haskey(data["param"], "dynamic_lambda")
            # Handle dynamic lambda
            gap = master.ext[:gap][end]
            threshold, EXTREME_GAP = data["param"]["dynamic_lambda"]
            old_lambda = master.ext[:stabilization_lambda][end]
            old_decr = master.ext[:stabilization_lambda_decr][end]

            if length(master.ext[:gap]) > 5 && all(master.ext[:gap][end-4:end] .> threshold * 0.5)
                lambda = min(1.0, old_lambda + decrement / 4)
                decr = decrement / 8
                phi = low_phi 
            elseif gap > 1.0
                # increment lambda back (retreat into stabilized region)
                lambda = min(1.0, old_lambda + old_decr / 2)
                # reduce the size of the decrementor
                decr = old_decr / 4
                phi = low_phi
            elseif gap > 0.02
                # keep lambda the same, explore and improve this region
                lambda = old_lambda
                decr = old_decr
                phi = low_phi
            else
                # decrease lambda (explore new areas more)
                lambda = max(0.0, old_lambda - old_decr)
                if lambda == 0 # if lambda would reach 0, just divide by 2 and converge to 0 instead
                    lambda = old_lambda / 2
                    decr = old_lambda / 2
                    phi = main_phi
                else
                    # increase the size of the decrementor
                    decr = old_decr * 2
                    phi = main_phi
                end
            end
        else
            # Handle linear lambda decrease
            lambda = max(min_lambda, starter_lambda - iter * decrement)
            phi = main_phi
        end
    end

    push!(master.ext[:stabilization_lambda], lambda)
    push!(master.ext[:stabilization_lambda_decr], decr)
    return lambda, phi
end

function compute_eval_core_points(master, data; trans_ratio=0.005, er_ratio=0.0005)
    iter = master.ext[:iter]

    trans_ratio = get(data["param"], "trans_perturb_ratio", trans_ratio)
    er_ratio = get(data["param"], "er_perturb_ratio", er_ratio)

    # Compute stabilization lambda
    stab_lambda, shift_phi = get_stabilization_shift(master, data)

    # core point perturbations
    gamma_eps = trans_ratio * data["param"]["num_cap_upgrades_max"]
    s_energy_eps = er_ratio * data["param"]["max_energy_rating"]

    y_core = nothing
    if iter == 0
        # Compute the core point
        gamma_core, s_energy_core = get_rep_day_core_point(data["param"]["core_point_simdir"])
        gamma_core = Float64.(gamma_core)
    
        gamma_core = max.(gamma_core, gamma_eps)
        s_energy_core = max.(s_energy_core, s_energy_eps)

        y_core = [gamma_core, s_energy_core]
        push!(master.ext[:y_core], y_core)

    else
        y_core_old = master.ext[:y_core][iter]
        y_eval_old = master.ext[:y_eval][iter]

        # Compute new shifted core point (convex combo of old core with old eval point)
        y_core = shift_phi * y_eval_old + (1 - shift_phi) * y_core_old
        push!(master.ext[:y_core], y_core)
    end

    # Compute the eval point (convex combo of raw with the core point)
    y_raw = master.ext[:y_raw][iter + 1]
    y_eval = stab_lambda * y_core + (1 - stab_lambda) * y_raw
    push!(master.ext[:y_eval], y_eval)

    return y_eval, y_core
end


function make_dual_subproblem(sub_model, data, y_eval)

    # Initialize model
    optimizer = Gurobi.Optimizer
    dual_sub_model = JuMP.Model(optimizer)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    # Tracked PTDF constraints in a dict [(a, r, t)] = true
    tracked_constraints = sub_model.ext[:tracked_constraints]
    # PTDF matrix
    PTDF = sub_model.ext[:PTDF]

    # Extract unique (arc, time) pairs from tracked constraints (ignoring representative r)
    arc_time_pairs = Set{Tuple{Int, Int}}()
    for (a, r, t) in keys(tracked_constraints)
        if sub_model.ext[:tracked_constraints][(a,r,t)]
            push!(arc_time_pairs, (a, t))
        end
    end

    # Unpack investment decisions
    gamma_eval, s_energy_eval = y_eval

    #
    #   I. Variables
    #

    # duals for power generation limits
    JuMP.@variable(dual_sub_model, mu_plus[g in 1:G, t in 1:T] <= 0)
    JuMP.@variable(dual_sub_model, mu_minus[g in 1:G, t in 1:T] >= 0)

    # duals for soc energy limits
    JuMP.@variable(dual_sub_model, tau_plus[i in 1:N, t in 1:T] <= 0)
    JuMP.@variable(dual_sub_model, tau_minus[i in 1:N, t in 1:T] >= 0)

    # duals for charging limits
    JuMP.@variable(dual_sub_model, psi_plus[i in 1:N, t in 1:T] <= 0)
    JuMP.@variable(dual_sub_model, psi_minus[i in 1:N, t in 1:T] >= 0)

    # dual for global power balance
    JuMP.@variable(dual_sub_model, alpha[t in 1:T])

    # dual for load-shedding non-negativity
    JuMP.@variable(dual_sub_model, beta[i in 1:N, t in 1:T] >= 0)

    # duals for initial/end soc
    JuMP.@variable(dual_sub_model, sigma_1[i in 1:N])
    JuMP.@variable(dual_sub_model, sigma_T[i in 1:N])

    # duals for soc tracking
    JuMP.@variable(dual_sub_model, rho[i in 1:N, t in 2:T])

    # duals for PTDF flow limits
    JuMP.@variable(dual_sub_model, pi_plus[(a, t) in arc_time_pairs] <= 0)
    JuMP.@variable(dual_sub_model, pi_minus[(a, t) in arc_time_pairs] >= 0)

    #
    #   II. Constraints
    #

    # dual for load_shedding ue
    JuMP.@constraint(dual_sub_model,
        ue[i=1:N, t=1:T],
        alpha[t] + beta[i, t] + sum(PTDF[a, i] * (pi_plus[(a, t)] + pi_minus[(a, t)]) 
                                    for (a, t_period) in arc_time_pairs if t_period == t) 
            == data["param"]["under_served_penalty"]
    )

    # dual for each generator's power
    JuMP.@constraint(dual_sub_model, 
        pg[g=1:G, t=1:T],
        mu_plus[g,t] + mu_minus[g,t] + alpha[t] + sum(PTDF[a, gen_bus_map[g]] * (pi_plus[(a, t)] + pi_minus[(a, t)]) 
                                        for (a, t_period) in arc_time_pairs if t_period == t)
            == data["gen"]["$g"]["cost"][2]
    )

    # dual for soc at hour 1
    JuMP.@constraint(dual_sub_model, 
        soc_1[i in 1:N],
        tau_plus[i,1] + tau_minus[i,1] - rho[i, 2] + sigma_1[i] == 0
    )

    # dual for soc at hours 2 - 23
    JuMP.@constraint(dual_sub_model,
        soc_middle[i in 1:N, t in 2:(T-1)],
        tau_plus[i,t] + tau_minus[i,t] + rho[i,t] - rho[i,t+1] == 0
    )

    # dual for soc at hour 24
    JuMP.@constraint(dual_sub_model,
        soc_T[i in 1:N],
        tau_plus[i,T] + tau_minus[i,T] + rho[i,T] + sigma_T[i] == 0
    )

    # dual for charge/discharge at hour 1
    JuMP.@constraint(dual_sub_model,
        ch_1[i in 1:N],
        psi_plus[i,1] + psi_minus[i,1] - alpha[1] - sigma_1[i] - sum(PTDF[a,i] * (pi_plus[(a,1)] + pi_minus[(a,1)]) 
            for (a, t_period) in arc_time_pairs if t_period == 1)
        == 0
    )

    # dual for charge/discharge at hour 2+
    JuMP.@constraint(dual_sub_model,
        ch[i in 1:N, t in 2:T],
        psi_plus[i,t] + psi_minus[i,t] - alpha[t] - rho[i,t] - sum(PTDF[a,i] * (pi_plus[(a,t)] + pi_minus[(a,t)])
            for (a, t_period) in arc_time_pairs if t_period == t)
        == 0
    )

    #
    #   III. Objective
    #
    @objective(dual_sub_model, Max, build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_eval, s_energy_eval))

    return dual_sub_model
end

function solve_for_pareto(dual_sub_model, sub_model, data, optimal_obj, y_eval, y_core)
    # Unpack the points
    gamma_eval, s_energy_eval = y_eval
    gamma_core, s_energy_core = y_core

    # Modify dual subproblem, add cut tightness constraint
    if haskey(data["param"], "cut_tightness") && data["param"]["cut_tightness"] == true
        eps = 1e-4 * abs(optimal_obj)
        @constraint(dual_sub_model,
            cut_tightness,
            build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_eval, s_energy_eval) >= optimal_obj - eps
        )
    end

    # Modify dual objective, to maximize violation at core point
    @objective(dual_sub_model, Max, build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_core, s_energy_core))
    optimize!(dual_sub_model)

    return dual_sub_model
end

function build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_vals, s_energy_vals)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    # Tracked PTDF constraints in a dict [(a, r, t)] = true
    tracked_constraints = sub_model.ext[:tracked_constraints]
    # PTDF matrix
    PTDF = sub_model.ext[:PTDF]

    # Extract unique (arc, time) pairs from tracked constraints (ignoring representative r)
    arc_time_pairs = Set{Tuple{Int, Int}}()
    for (a, r, t) in keys(tracked_constraints)
        if sub_model.ext[:tracked_constraints][(a,r,t)]
            push!(arc_time_pairs, (a, t))
        end
    end

    # Access variables from the model
    mu_plus = dual_sub_model[:mu_plus]
    mu_minus = dual_sub_model[:mu_minus]
    tau_plus = dual_sub_model[:tau_plus]
    psi_plus = dual_sub_model[:psi_plus]
    psi_minus = dual_sub_model[:psi_minus]
    alpha = dual_sub_model[:alpha]
    pi_plus = dual_sub_model[:pi_plus]
    pi_minus = dual_sub_model[:pi_minus]
    sigma_1 = dual_sub_model[:sigma_1]
    sigma_T = dual_sub_model[:sigma_T]

    # Objective function
    operational_weight = get(data["param"], "operational_weight", 1)

    pg_var = sub_model[:pg]
    pg_upper_bounds = upper_bound.(pg_var)
    pg_lower_bounds = lower_bound.(pg_var)

    return operational_weight * 
        (sum(sum(
                mu_plus[g,t] * pg_upper_bounds[1, g, t]
                + mu_minus[g,t] * pg_lower_bounds[1, g, t]
            for t in 1:T)
        for g in 1:G)
        + sum(sum(
                tau_plus[i,t] * s_energy_vals[i]
                + psi_plus[i,t] * s_power_vals[i]
                - psi_minus[i,t] * s_power_vals[i]
                + alpha[t] * data["bus"]["$i"]["load"]["1"][t]
            for i in 1:N)
        for t in 1:T)
        + sum(
            get(data["param"], "soc_init_end_ratio", 0.5) * s_energy_vals[i] * (sigma_1[i] + sigma_T[i])
        for i in 1:N)
        + sum(
            pi_plus[(a,t)] * (
                data["branch"]["$a"]["rate_a"]
                + get_capacity_increment(data, a) * gamma_vals[a]
                + sum(
                    PTDF[a,i] * data["bus"]["$i"]["load"]["1"][t]
                for i in 1:N)
            )
            + pi_minus[(a,t)] * (
                - data["branch"]["$a"]["rate_a"]
                - get_capacity_increment(data, a) * gamma_vals[a]
                + sum(
                    PTDF[a,i] * data["bus"]["$i"]["load"]["1"][t]
                for i in 1:N)
            )
        for (a,t) in arc_time_pairs))
end

function add_pareto_cut(master, dual_sub_model, sub_model, data, theta, y)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    # Tracked PTDF constraints in a dict [(a, r, t)] = true
    tracked_constraints = sub_model.ext[:tracked_constraints]
    # PTDF matrix
    PTDF = sub_model.ext[:PTDF]

    # Extract unique (arc, time) pairs from tracked constraints (ignoring representative r)
    arc_time_pairs = Set{Tuple{Int, Int}}()
    for (a, r, t) in keys(tracked_constraints)
        if sub_model.ext[:tracked_constraints][(a,r,t)]
            push!(arc_time_pairs, (a, t))
        end
    end

    # Unpack y variable
    gamma, s_power, s_energy = y

    # Get other model parameters
    operational_weight = get(data["param"], "operational_weight", 1)
    pg_var = sub_model[:pg]
    pg_upper_bounds = upper_bound.(pg_var)
    pg_lower_bounds = lower_bound.(pg_var)

    # Unpack the optimal duals from the solved dual_sub_model
    mu_plus = value.(dual_sub_model[:mu_plus])
    mu_minus = value.(dual_sub_model[:mu_minus])
    tau_plus = value.(dual_sub_model[:tau_plus])
    psi_plus = value.(dual_sub_model[:psi_plus])
    psi_minus = value.(dual_sub_model[:psi_minus])
    alpha = value.(dual_sub_model[:alpha])
    sigma_1 = value.(dual_sub_model[:sigma_1])
    sigma_T = value.(dual_sub_model[:sigma_T])
    pi_plus = value.(dual_sub_model[:pi_plus])
    pi_minus = value.(dual_sub_model[:pi_minus])

    @constraint(master,
        theta >= operational_weight * 
        (sum(sum(
                mu_plus[g,t] * pg_upper_bounds[1, g, t]
                + mu_minus[g,t] * pg_lower_bounds[1, g, t]
            for t in 1:T)
        for g in 1:G)
        + sum(sum(
                tau_plus[i,t] * s_energy[i]
                + psi_plus[i,t] * s_power[i]
                - psi_minus[i,t] * s_power[i]
                + alpha[t] * data["bus"]["$i"]["load"]["1"][t]
            for i in 1:N)
        for t in 1:T)
        + sum(
            get(data["param"], "soc_init_end_ratio", 0.5) * s_energy[i] * (sigma_1[i] + sigma_T[i])
        for i in 1:N)
        + sum(
            pi_plus[(a,t)] * (
                data["branch"]["$a"]["rate_a"]
                + get_capacity_increment(data, a) * gamma[a]
                + sum(
                    PTDF[a,i] * data["bus"]["$i"]["load"]["1"][t]
                for i in 1:N)
            )
            + pi_minus[(a,t)] * (
                - data["branch"]["$a"]["rate_a"]
                - get_capacity_increment(data, a) * gamma[a]
                + sum(
                    PTDF[a,i] * data["bus"]["$i"]["load"]["1"][t]
                for i in 1:N)
            )
        for (a,t) in arc_time_pairs))
    )
end


    