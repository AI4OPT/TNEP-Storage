function solve_subrel_transport(simdir, y_val, data)
    # unpack investment decisions
    gamma_val, s_power_val, s_energy_val = y_val

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])

    # Initialize model
    subrel_optimizer = Gurobi.Optimizer
    subrel = JuMP.Model(subrel_optimizer)
    set_optimizer_attribute(subrel, "LogFile", joinpath(simdir, "gurobi_subrel_logfile.log"))
    set_optimizer_attribute(subrel, "MIPGap", data["param"]["mip_gap"])

    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    
    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))
    
    # Branch metadata
    from_bus = Dict(parse(Int, a) => data["branch"]["$a"]["f_bus"] for a in keys(data["branch"]))
    to_bus = Dict(parse(Int, a) => data["branch"]["$a"]["t_bus"] for a in keys(data["branch"]))
    
    #
    #   I. Variables
    #
    
    # Generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    JuMP.@variable(subrel, pg[r in 1:R, g in 1:G, t in 1:T])

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

    # Under-served energy
    @variable(subrel, ue[r=1:R, i=1:N, t=1:T] >= 0)
    
    # Storage operation
    @variable(subrel, ch[r=1:R, i=1:N, t=1:T])  # Charging/discharging
    @variable(subrel, soc[r=1:R, i=1:N, t=1:T] >= 0)  # State of charge

    # Power flow variables (transport model)
    @variable(subrel, flow[r=1:R, a=1:E, t=1:T])
    
    # Fix investment variables based on master problem decisions
    @variable(subrel, gamma[a=1:E])
    @variable(subrel, s_power[i=1:N])
    @variable(subrel, s_energy[i=1:N])

    #
    #   II. Constraints
    #

    # FIX INVESTMENT DECISIONS FROM MASTER PROBLEM (DUALS USED FOR CUTS)
    @constraint(subrel, master_gamma[a=1:E], gamma[a] == gamma_val[a])
    @constraint(subrel, master_power[i=1:N], s_power[i] == s_power_val[i])
    @constraint(subrel, master_energy[i=1:N], s_energy[i] == s_energy_val[i])

    # Node power balance (conservation of flow)
    @constraint(subrel, 
        node_balance[r=1:R, i=1:N, t=1:T],
        sum(pg[r,g,t] for g=1:G if gen_bus_map[g] == i; init=0.0) + 
        sum(flow[r,a,t] for a=1:E if to_bus[a] == i; init=0.0) -
        sum(flow[r,a,t] for a=1:E if from_bus[a] == i; init=0.0) -
        data["bus"]["$i"]["load"]["$r"][t] + 
        ue[r,i,t] - 
        ch[r,i,t] == 0
    )

    # Line capacity constraints
    @constraint(subrel,
        line_capacity_ub[r=1:R, a=1:E, t=1:T],
        flow[r,a,t] <= data["branch"]["$a"]["rate_a"] + gamma[a] * get_capacity_increment(data, a)
    )
    
    @constraint(subrel,
        line_capacity_lb[r=1:R, a=1:E, t=1:T],
        flow[r,a,t] >= -(data["branch"]["$a"]["rate_a"] + gamma[a] * get_capacity_increment(data, a))
    )

    # Storage constraints

    # SOC evolution
    @constraint(subrel, 
        soc_over_time[r=1:R, i=1:N, t=2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t]
    )
    
    # Initial and final SOC
    @constraint(subrel,
        soc_start[r=1:R, i=1:N],
        soc[r,i,1] == 0.5 * s_energy[i] + ch[r,i,1]
    )
    @constraint(subrel,
        soc_end[r=1:R, i=1:N],
        soc[r,i,T] == 0.5 * s_energy[i]
    )

    # SOC limits
    @constraint(subrel, 
        soc_energy_ub[r=1:R, i=1:N, t=1:T],
        soc[r,i,t] <= s_energy[i]
    )
    
    # Charging/discharging limits
    @constraint(subrel,
        charge_discharge_lb[r=1:R, i=1:N, t=1:T],
        -s_power[i] <= ch[r,i,t]
    )
    @constraint(subrel,
        charge_discharge_ub[r=1:R, i=1:N, t=1:T],
        ch[r,i,t] <= s_power[i]
    )

    # Objective function
    operational_weight = get(data["param"], "operational_weight", 1)

    @objective(subrel, Min,
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

    # Solve the model
    t0 = time()
    optimize!(subrel)
    solve_time = time() - t0
    
    # Check solution status
    if termination_status(subrel) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        error("Transport subrelproblem failed to solve optimally: $(termination_status(subrel))")
    end
    
    # Record solve time
    subrel.ext[:solve_time] = solve_time
    println("Transport subrelproblem solved in $solve_time seconds")
    
    # Extract duals for Benders cuts
    dual_gamma = [dual(master_gamma[a]) for a=1:E]
    dual_power = [dual(master_power[i]) for i=1:N]
    
    duals = (dual_gamma, dual_power)
    
    # Return solution
    return subrel, objective_value(subrel), duals
end