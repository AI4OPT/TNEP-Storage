using JuMP
using Gurobi
include("../../helpers/compute_gen_cost.jl")
include("naive_candidates.jl")
include("rate_a_zero.jl")

function define_master_transmission(data::Dict{String, Any}; prev_simdir=nothing)

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

    # investment level of capacity upgrade
    JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)

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

    #
    #   III. Objective
    #

    JuMP.@objective(master, Min,
    sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
    theta
    )

    y = gamma
    return master, y, theta
end

function solve_subproblem_transmission(simdir, y_val, data)

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

    # investment level of capacity upgrade
    JuMP.@variable(sub, gamma[a=1:E])

    # branch flows
    JuMP.@variable(sub, pf[r=1:R, a=1:E, t=1:T])

    # voltage angles
    JuMP.@variable(sub, va[r=1:R, i=1:N, t=1:T])

    #
    #   II. Constraints
    #

    # flow constraints
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    JuMP.@constraint(sub, 
        flow_lb[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] >= -1 * data["branch"]["$a"]["rate_a"] - gamma[a] * data["param"]["cap_upgrade_increment"]
    )
    JuMP.@constraint(sub, 
        flow_ub[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] <= data["branch"]["$a"]["rate_a"] + gamma[a] * data["param"]["cap_upgrade_increment"]
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
        ue[r,i,t] ==
        sum(pf[r,a,t] for a in [a for a in data["arcs_from"]["$i"] if data["branch"]["$a"]["f_bus"] == i]) + # outflow
        data["bus"]["$i"]["load"]["$r"][t] +
        oe[r,i,t]
    )

    # MASTER PROBLEM FIXED DECISIONS
    gamma_val = y_val
    JuMP.@constraint(sub,
        master_gamma[a in 1:E],
        gamma[a] == gamma_val[a]
    )

    #
    #   III. Objective
    #
    operational_weight = 1
    if haskey(data["param"], "operational_weight")
        operational_weight = data["param"]["operational_weight"]
    end

    JuMP.@objective(sub, Min,
    sum(
        data["param"]["representative_prob"]["$r"] *
        (
            sum(
                sum(compute_gen_cost(pg[r,g,t], data["gen"]["$g"]) for g in 1:G) +
                sum(data["param"]["over_generated_penalty"] * oe[r,i,t] for i in 1:N) + 
                sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
            for t in 1:T)
        )
    for r in 1:R) * operational_weight
    )

    # Optimize and return duals
    set_optimizer_attribute(sub, "LogFile", joinpath(simdir, "gurobi_sub_logfile.log"))
    set_optimizer_attribute(sub, "MIPGap", data["param"]["mip_gap"])
    optimize!(sub)

    dual_master_gamma = [dual(master_gamma[a]) for a in 1:E]
    duals = dual_master_gamma
    return objective_value(sub), duals
end

function add_benders_cut_transmission(master, theta, duals, y, y_val, phi_val)
    JuMP.@constraint(master, 
    theta >= phi_val + sum(duals .* (y .- y_val)))
end

function benders_iteration_transmission(simdir, master, y, theta, data)
    converged = false
    iter = 0
    gamma_val = nothing

    while !converged
        # solve master and decide investments
        optimize!(master)
        gamma_val = value.(master[:gamma])
        y_val = gamma_val
        theta_val = value.(master[:theta])
        
        # get the subproblem objective and duals
        phi_val, duals = solve_subproblem_transmission(simdir, y_val, data)

        # save the progress in a csv
        filename = joinpath(simdir, "output", "benders_progress.csv")
        benders_transmission_write_to_csv(filename, objective_value(master), theta_val, phi_val, gamma_val)

        # check for convergence, otherwise, add a new cut to the master problem
        if is_converged(phi_val, theta_val)
            converged = true
        else
            y = master[:gamma]
            theta = master[:theta]
            add_benders_cut_transmission(master, theta, duals, y, y_val, phi_val)
        end

        iter += 1
    end

    return gamma_val
end

function benders_only_transmission(simdir, data)
    master, y, theta = define_master_transmission(data)
    gamma_val = benders_iteration_transmission(simdir, master, y, theta, data)

    df = DataFrame(gamma = gamma_val)
    filename = joinpath(simdir, "output", "investments.csv")
    CSV.write(filename, df)
    return
end

function is_converged(phi_val, theta_val, tol=0.01)
    if phi_val == 0
        return abs(theta_val) < tol * (1 + abs(theta_val))
    else
        return abs(theta_val - phi_val) / abs(phi_val) < tol
    end
end

function benders_transmission_write_to_csv(filename, master_obj, theta_val, phi_val, gamma_val)
    # Prepare the data row
    df = DataFrame(
        master_objective = [master_obj],
        theta_val = [theta_val],
        phi_val = [phi_val]
    )

    # Write to CSV, appending to it
    if isfile(filename)
        CSV.write(filename, df; append = true, header = false)
    else
        CSV.write(filename, df)
    end
end



