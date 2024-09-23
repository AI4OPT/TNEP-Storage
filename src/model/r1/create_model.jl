# input data will be PowerModels dictionary
# slightly modified to contain more information

"""
data {
    param: {
        num_cap_upgrades_max:
        num_representatives:
        num_hours:
        cap_upgrade_cost:
        ...
    }

    scenario_prob: {
        1: 0.50
        2: 0.50
    }

    gen: {
       1: {
        gen_type: wind
        ...
        profile: {
            1: [100, 90, 110, 120]
            2: [80, 90, 70, 100]
        }
       }
       2: {
        gen_type: thermal
        ...
        }
       }
       ...
    }

    bus: {
        1: {
            bus_name: 
            ...
            load: {
                1: [60, 70, 100, 90]
            },
            gen: {thermal: [1], wind: [2], solar: [3]}
        }
    },
    arcs_from: {                # key: bus_index => value: branch_index
        "1": ["3", "5"]
    },

    branch: {
        "3": {
            f_bus: 1
            t_bus: 7
            ...
            "cap_upgrade_increment":
            "cap_upgrade_cost":
        },

        "5": {
            f_bus: 9
            t_bus: 1
            ...
        }
    }
}
"""

using JuMP
include("../../helpers/compute_gen_cost.jl")
include("naive_candidates.jl")
include("rate_a_zero.jl")

function create_model_r1(data::Dict{String, Any}, optimizer; prev_simdir=nothing)

    # Initialize model
    model = JuMP.Model(optimizer)

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
    @variable(model, pg[r in 1:R, g in 1:G, t in 1:T])

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
    JuMP.@variable(model, oe[r=1:R, i=1:N, t=1:T] >= 0)

    # under-served energy at bus
    JuMP.@variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0)

    # investment level of capacity upgrade
    JuMP.@variable(model, 0 <= gamma[a=1:E] <= K, Int)

    # branch flows
    JuMP.@variable(model, pf[r=1:R, a=1:E, t=1:T])

    # voltage angles
    JuMP.@variable(model, va[r=1:R, i=1:N, t=1:T])

    # power rating of storage
    JuMP.@variable(model, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(model, s_energy[i=1:N] >= 0)

    # state of charge of storage
    JuMP.@variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0)

    # charging of storage
    JuMP.@variable(model, ch[r=1:R, i=1:N, t=1:T] >= 0)

    # discharging of storage
    JuMP.@variable(model, dis[r=1:R, i=1:N, t=1:T] >= 0)

    # binary variable for installation of storage
    JuMP.@variable(model, sigma[i=1:N], Bin)

    # binary variable indicating charge/discharge
    JuMP.@variable(model, alpha[r=1:R, i=1:N, t=1:T], Bin)
    
    # linearizing (continuous) variable
    JuMP.@variable(model, beta[r=1:R, i=1:N, t=1:T] >= 0)

    #
    #   II. Constraints
    #

    # if prev model as input, enforce old investment upgrades (gamma, s_power, s_energy)
    if prev_simdir !== nothing
        line_inv = CSV.read(joinpath(prev_simdir, "output", "line_investments.csv"), DataFrame)
        storage_inv = CSV.read(joinpath(prev_simdir, "output", "storage_investments.csv"), DataFrame)
        JuMP.@constraint(model,
            old_gamma[a in 1:E],
            gamma[a] >= line_inv[a, :Upgrade_Lvl]
        )
        JuMP.@constraint(model,
            old_s_power[i in 1:N],
            s_power[i] >= storage_inv[i, :Storage_Power]
        )
        JuMP.@constraint(model,
            old_s_energy[i in 1:N],
            s_energy[i] >= storage_inv[i, :Storage_Energy]
        )
    end

    # flow constraints
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    JuMP.@constraint(model, 
        flow_lb[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] >= -1 * data["branch"]["$a"]["rate_a"] - gamma[a] * data["param"]["cap_upgrade_increment"]
    )
    JuMP.@constraint(model, 
        flow_ub[r in 1:R, a in rate_a_nonzero, t in 1:T],
        pf[r,a,t] <= data["branch"]["$a"]["rate_a"] + gamma[a] * data["param"]["cap_upgrade_increment"]
    )

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # Ohm's law constraint
    JuMP.@constraint(model,
        ohms_law[r in 1:R, a in 1:E, t in 1:T],
        pf[r,a,t] == (va[r,tbus[a],t] - va[r,fbus[a],t]) / br_x[a]
    )

    # Optional: slack bus voltage angle
    JuMP.@constraint(model,
    slack_bus_voltage[r in 1:R, t in 1:T],
    va[r, 1, t] == 0
    )

    # Voltage difference limits
    JuMP.@constraint(model,
        voltage_diff_ub[r in 1:R, a in 1:E, t in 1:T],
        va[r,tbus[a],t] - va[r,fbus[a],t] <= data["param"]["voltage_angle_difference_max"]
    )
    JuMP.@constraint(model,
        voltage_diff_lb[r in 1:R, a in 1:E, t in 1:T],
        va[r,tbus[a],t] - va[r,fbus[a],t] >= -1 * data["param"]["voltage_angle_difference_max"]
    )

    # nodal power balance
    JuMP.@constraint(model, 
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
    JuMP.@constraint(model, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t] * data["param"]["bess_efficiency"] - dis[r,i,t] / data["param"]["bess_efficiency"]
    )

    # OPTIONAL: soc 0.5 constraint
    JuMP.@constraint(model,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == 0.5 * s_energy[i] + ch[r,i,1] * data["param"]["bess_efficiency"] - dis[r,i,1] / data["param"]["bess_efficiency"]
    )
    JuMP.@constraint(model,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == 0.5 * s_energy[i]
    )

    # OPTIONAL: limit the number of storage locations
    if (haskey(data["param"], "max_num_storage") && data["param"]["max_num_storage"] >= 0)
        JuMP.@constraint(model, 
            storage_location_limit,
            sum(sigma[i] for i in 1:N) <= data["param"]["max_num_storage"]
        )
    end

    # soc energy rating constraint
    JuMP.@constraint(model, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i]
    )

    # energy rating only if storage installed
    JuMP.@constraint(model, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= sigma[i] * data["param"]["max_energy_rating"]
    )

    # power rating only if storage installed
    JuMP.@constraint(model, 
        installed_power_ub[i in 1:N],
        s_power[i] <= sigma[i] * data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    JuMP.@constraint(model, 
        short_duration[i in 1:N],
        s_energy[i] <= 4.0 * s_power[i]
    )

    # ~~~ STORAGE SOC LINEARIZATION CONSTRAINTS ~~~
    JuMP.@constraint(model, 
        storage_linearization_1[r in 1:R, i in 1:N, t in 1:T],
        data["param"]["bess_efficiency"] * ch[r,i,t] <= beta[r,i,t]
    )
    JuMP.@constraint(model, 
        storage_linearization_2[r in 1:R, i in 1:N, t in 1:T],
        dis[r,i,t] / data["param"]["bess_efficiency"] <= s_power[i] - beta[r,i,t]
    )
    JuMP.@constraint(model, 
        storage_linearization_3[r in 1:R, i in 1:N, t in 1:T],
        beta[r,i,t] <= data["param"]["max_power_rating"] * alpha[r,i,t]
    )
    JuMP.@constraint(model, 
        storage_linearization_4[r in 1:R, i in 1:N, t in 1:T],
        s_power[i] - beta[r,i,t] <= data["param"]["max_power_rating"] * (1 - alpha[r,i,t])
    )

    # OPTIONAL: CANDIDATE STORAGE LOCATIONS ONLY
    if haskey(data["param"], "candidate_analysis_seqsimdir")
        seqsimdir = data["param"]["candidate_analysis_seqsimdir"]
        decarbonization_year = data["param"]["decarbonization_year"]
        filepath = "$(seqsimdir)/$(decarbonization_year)/output/energy.csv"
        candidates, non_candidates = get_storage_candidates(data, filepath)
        non_candidates = Set(parse(Int, x) for x in non_candidates)

        # JuMP.@constraint(model, storage_non_candidate[i in non_candidates],
        #     sigma[i] == 0
        # )
        for i in non_candidates
            fix(sigma[i], 0; force = true)
        end
    end

    #
    #   III. Objective
    #

    operational_weight = 1
    if haskey(data["param"], "operational_weight")
        operational_weight = data["param"]["operational_weight"]
    end

    JuMP.@objective(model, Min,
    sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
    sum(sigma[i] for i in 1:N) * data["param"]["storage_fixed_cost"] + 
    sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
    sum(
        data["param"]["representative_prob"][r] *
        (
            sum(
                sum(compute_gen_cost(pg[r,g,t], data["gen"]["$g"]) for g in 1:G) +
                sum(data["param"]["over_generated_penalty"] * oe[r,i,t] for i in 1:N) + 
                sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
            for t in 1:T)
        )
    for r in 1:R) * operational_weight
    )
    
    return model
end