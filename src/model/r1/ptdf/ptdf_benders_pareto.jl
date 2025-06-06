function make_dual_subproblem(sub_model, data, y_val)

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
    gamma_val, s_power_val, s_energy_val = y_val

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
    @objective(dual_sub_model, Max, build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_val, s_power_val, s_energy_val))

    return dual_sub_model
end

function solve_for_pareto(dual_sub_model, data, optimal_obj, master_y_val, core_y_val)
    operational_weight = get(data["param"], "operational_weight", 1)

    # Unpack the points
    master_gamma_val, master_s_power_val, master_s_energy_val = master_y_val
    cor_gamma_val, cor_s_power_val, cor_s_energy_val = core_y_val

    # Modify dual subproblem, add cut validity constraint
    @constraint(dual_sub_model,
        cut_validity,
        build_dual_obj_expr(dual_sub_model, sub_model, data, master_gamma_val, master_s_power_val, master_s_energy_val) == optimal_obj
    )

    # Modify dual objective, to maximize violation at core point
    @objective(dual_sub_model, Max, build_dual_obj_expr(dual_sub_model, sub_model, data, cor_gamma_val, cor_s_power_val, cor_s_energy_val))
    optimize!(dual_sub_model)

    return dual_sub_model
end

function build_dual_obj_expr(dual_sub_model, sub_model, data, gamma_vals, s_power_vals, s_energy_vals)
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
            0.5 * s_energy_vals[i] * (sigma_1[i] + sigma_T[i])
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
            0.5 * s_energy[i] * (sigma_1[i] + sigma_T[i])
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
    