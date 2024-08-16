using JuMP
using Gurobi
include("../../helpers/compute_gen_cost.jl")
include("naive_candidates.jl")
include("rate_a_zero.jl")

function define_master_storage(data::Dict{String, Any}; prev_simdir=nothing)

    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)

    # Establish index sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    fbus::Vector{Int} = [data["branch"]["$a"]["f_bus"] for a in 1:E]                # from bus
    tbus::Vector{Int} = [data["branch"]["$a"]["t_bus"] for a in 1:E]                # to bus
    br_x::Vector{Float64} = [data["branch"]["$a"]["br_x"] for a in 1:E]             # branch reactance

    #
    #   I. Variables
    #

    # power rating of storage
    JuMP.@variable(master, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(master, s_energy[i=1:N] >= 0)

    # binary variable for installation of storage
    JuMP.@variable(master, sigma[i=1:N], Bin)

    # subproblem objective(s)
    JuMP.@variable(master, theta >= 0)

    #
    #   II. Constraints
    #

    # energy rating only if storage installed
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= sigma[i] * data["param"]["max_energy_rating"]
    )

    # power rating only if storage installed
    JuMP.@constraint(master, 
        installed_power_ub[i in 1:N],
        s_power[i] <= sigma[i] * data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    JuMP.@constraint(master, 
        short_duration[i in 1:N],
        s_energy[i] <= 4.0 * s_power[i]
    )

    # OPTIONAL: CANDIDATE STORAGE LOCATIONS ONLY
    if haskey(data["param"], "candidate_analysis_seqsimdir")
        seqsimdir = data["param"]["candidate_analysis_seqsimdir"]
        decarbonization_year = data["param"]["decarbonization_year"]
        filepath = "$(seqsimdir)/$(decarbonization_year)/output/energy.csv"
        candidates, non_candidates = get_storage_candidates(data, filepath)
        non_candidates = Set(parse(Int, x) for x in non_candidates)

        # JuMP.@constraint(master, storage_non_candidate[i in non_candidates],
        #     sigma[i] == 0
        # )
        for i in non_candidates
            fix(sigma[i], 0; force = true)
        end
    end

    #
    #   III. Objective
    #
    JuMP.@objective(master, Min,
    sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
    sum(sigma[i] for i in 1:N) * data["param"]["storage_fixed_cost"] +
    theta
    )

    y = s_power, s_energy, sigma

    return master, y, theta
end

function solve_subproblem_storage(simdir, y_val, data)

    # Initialize model
    sub_optimizer = Gurobi.Optimizer
    sub = JuMP.Model(sub_optimizer)

    # Establish index sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    fbus::Vector{Int} = [data["branch"]["$a"]["f_bus"] for a in 1:E]                # from bus
    tbus::Vector{Int} = [data["branch"]["$a"]["t_bus"] for a in 1:E]                # to bus
    br_x::Vector{Float64} = [data["branch"]["$a"]["br_x"] for a in 1:E]             # branch reactance

    #
    #   I. Variables
    #

    # generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    @variable(sub, pg[r in 1:R, g in 1:G, t in 1:T])

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

    # over-generated energy at bus
    JuMP.@variable(sub, oe[r=1:R, i=1:N, t=1:T] >= 0)

    # under-served energy at bus
    JuMP.@variable(sub, ue[r=1:R, i=1:N, t=1:T] >= 0)

    # branch flows
    JuMP.@variable(sub, pf[r=1:R, a=1:E, t=1:T])

    # power rating of storage
    JuMP.@variable(sub, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(sub, s_energy[i=1:N] >= 0)

    # voltage angles
    JuMP.@variable(sub, va[r=1:R, i=1:N, t=1:T])

    # state of charge of storage
    JuMP.@variable(sub, soc[r=1:R, i=1:N, t=1:T] >= 0)

    # charging of storage
    JuMP.@variable(sub, ch[r=1:R, i=1:N, t=1:T] >= 0)

    # discharging of storage
    JuMP.@variable(sub, dis[r=1:R, i=1:N, t=1:T] >= 0)

    # binary variable for installation of storage
    JuMP.@variable(sub, sigma[i=1:N], Bin)

    # binary variable indicating charge/discharge
    JuMP.@variable(sub, alpha[r=1:R, i=1:N, t=1:T], Bin)
    
    # linearizing (continuous) variable
    JuMP.@variable(sub, beta[r=1:R, i=1:N, t=1:T] >= 0)

    #
    #   II. Constraints
    #

    # flow constraints
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    JuMP.@constraint(sub, 
        flow_lb[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] >= -1 * data["branch"]["$a"]["rate_a"]
    )
    JuMP.@constraint(sub, 
        flow_ub[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] <= data["branch"]["$a"]["rate_a"]
    )

    # Ohm's law constraint
    JuMP.@constraint(sub,
        ohms_law[r in 1:R, a in 1:E, t in 1:T],
        pf[r,a,t] == (va[r,tbus[a],t] - va[r,fbus[a],t]) / br_x[a]
    )

    # Optional: slack bus voltage angle
    JuMP.@constraint(sub,
    slack_bus_voltage[r in 1:R, t in 1:T],
    va[r, 1, t] == 0
    )

    # Voltage difference limits
    JuMP.@constraint(sub,
        voltage_diff_ub[r in 1:R, a in 1:E, t in 1:T],
        va[r,tbus[a],t] - va[r,fbus[a],t] <= data["param"]["voltage_angle_difference_max"]
    )
    JuMP.@constraint(sub,
        voltage_diff_lb[r in 1:R, a in 1:E, t in 1:T],
        va[r,tbus[a],t] - va[r,fbus[a],t] >= -1 * data["param"]["voltage_angle_difference_max"]
    )

    # nodal power balance
    JuMP.@constraint(sub, 
        power_balance[r in 1:R, i in 1:N, t in 1:T],
        sum(pg[r,g,t] for g in [num for gen_type in values(data["bus"]["$i"]["gen"]) if isa(gen_type, Array) for num in gen_type]) +
        sum(pf[r,a,t] for a in [a for a in data["arcs_from"]["$i"] if data["branch"]["$a"]["t_bus"] == i]) + # inflow
        dis[r,i,t] +
        ue[r,i,t] ==
        sum(pf[r,a,t] for a in [a for a in data["arcs_from"]["$i"] if data["branch"]["$a"]["f_bus"] == i]) + # outflow
        data["bus"]["$i"]["load"]["$r"][t] +
        oe[r,i,t] +
        ch[r,i,t]
    )

    # soc over time constraint
    JuMP.@constraint(sub, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - dis[r,i,t] / data["param"]["bess_efficiency"]
    )

    # OPTIONAL: soc 0.5 constraint
    JuMP.@constraint(sub,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == 0.5 * s_energy[i] + ch[r,i,1] * data["param"]["bess_efficiency"] - dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    JuMP.@constraint(sub,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == 0.5 * s_energy[i]
    )

    # soc energy rating constraint
    JuMP.@constraint(sub, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i]
    )

    # ~~~ STORAGE SOC LINEARIZATION CONSTRAINTS ~~~
    JuMP.@constraint(sub, 
        storage_linearization_1[r in 1:R, i in 1:N, t in 1:T],
        data["param"]["bess_efficiency"] * ch[r,i,t] <= beta[r,i,t]
    )
    JuMP.@constraint(sub, 
        storage_linearization_2[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] / data["param"]["bess_efficiency"] <= s_power[i] - beta[r,i,t]
    )
    JuMP.@constraint(sub, 
        storage_linearization_3[r in 1:R, i in 1:N, t in 1:T],
        beta[r,i,t] <= data["param"]["max_power_rating"] * alpha[r,i,t]
    )
    JuMP.@constraint(sub, 
        storage_linearization_4[r in 1:R, i in 1:N, t in 1:T],
        s_power[i] - beta[r,i,t] <= data["param"]["max_power_rating"] * (1 - alpha[r,i,t])
    )

    # MASTER PROBLEM FIXED DECISIONS
    s_power_val, s_energy_val, sigma_val = y_val
    JuMP.@constraint(sub,
        master_power[i in 1:N],
        s_power[i] == s_power_val[i]
    )
    JuMP.@constraint(sub,
        master_energy[i in 1:N],
        s_energy[i] == s_energy_val[i]
    )
    JuMP.@constraint(sub,
        master_sigma[i in 1:N],
        sigma[i] == sigma_val[i]
    )

    #
    #   III. Objective
    #

    JuMP.@objective(sub, Min,
    sum(
        data["representative_prob"]["$r"] *
        (
            sum(
                sum(compute_gen_cost(pg[r,g,t], data["gen"]["$g"]) for g in 1:G) +
                sum(data["param"]["over_generated_penalty"] * oe[r,i,t] for i in 1:N) + 
                sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
            for t in 1:T)
        )
    for r in 1:R)
    )

    # Optimize and return duals
    set_optimizer_attribute(sub, "LogFile", joinpath(simdir, "gurobi_sub_logfile.log"))
    set_optimizer_attribute(sub, "MIPGap", data["param"]["mip_gap"])
    optimize!(sub)

    dual_master_power = [dual(master_power[i]) for i in 1:N]
    dual_master_energy = [dual(master_energy[i]) for i in 1:N]
    dual_master_sigma = [dual(master_sigma[i]) for i in 1:N]
    duals = dual_master_power, dual_master_energy, dual_master_sigma
    return objective_value(sub), duals
end

function add_benders_cut_storage(master, theta, duals, y, y_val, theta_val)
    y_concat = vcat(y[1], y[2], y[3])   # Combine s_power, s_energy, sigma
    y_val_concat = vcat(y_val[1], y_val[2], y_val[3])
    JuMP.@constraint(master, theta >= theta_val + dot(duals, y .- y_val))
end

function benders_iteration_storage(simdir, master, y, theta, data)
    converged = false
    while !converged
        # solve master and decide investments
        optimize!(master)
        s_power_val = value.(master[:s_power])
        s_energy_val = value.(master[:s_energy])
        sigma_val = value.(master[:sigma])
        y_val = s_power_val, s_energy_val, sigma_val
        
        # get the subproblem objective and duals
        theta_val, duals = solve_subproblem_storage(simdir, y_val, data)
        
        if is_converged(theta_val, value.(master[:theta]))
            converged = true
        else
            y = s_power, s_energy, sigma
            add_benders_cut_storage(master, theta, duals, y, y_val, theta_val)
        end
    end
end

function benders_only_storage(simdir, data)
    master, y, theta = define_master_storage(data)
    benders_iteration_storage(simdir, master, y, theta, data)