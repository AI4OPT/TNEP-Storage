using JuMP
include("../../../helpers/compute_gen_cost.jl")
include("../storage_candidates/naive_candidates.jl")
include("../rate_a_zero.jl")
include("base_ptdf.jl")

function create_model_r1_ptdf(data::Dict{String, Any}, optimizer; line_investments=nothing, storage_investments=nothing)

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

    # PTDF matrix
    ptdf = do_all_ptdf(data)

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

    # under-served energy at bus
    JuMP.@variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0)

    # investment level of capacity upgrade
    JuMP.@variable(model, 0 <= gamma[a=1:E] <= K, Int)

    # branch flows
    JuMP.@variable(model, pf[r=1:R, a=1:E, t=1:T])

    # power rating of storage
    JuMP.@variable(model, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(model, s_energy[i=1:N] >= 0)

    # state of charge of storage
    JuMP.@variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0)

    # charging of storage
    JuMP.@variable(model, ch[r=1:R, i=1:N, t=1:T])

    # binary variable for installation of storage
    JuMP.@variable(model, sigma[i=1:N], Bin)

    #
    #   II. Constraints
    #

    # if prev model as input, enforce old investment upgrades (gamma, s_power, s_energy)
    nonzero_storage_nodes = Set()
    if haskey(data["param"], "prev_simdir")
        prev_simdir = data["param"]["prev_simdir"]
        line_inv = CSV.read(joinpath(prev_simdir, "output", "line_investments.csv"), DataFrame)
        storage_inv = CSV.read(joinpath(prev_simdir, "output", "storage_investments.csv"), DataFrame)
        for i in 1:N
            if storage_inv[i, :Storage_Power] != 0 || storage_inv[i, :Storage_Energy] != 0
                push!(nonzero_storage_nodes, "$i")
            end
        end
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

    # if there are investments to test, enforce these exact upgrades
    if line_investments !== nothing
        inv = CSV.read(line_investments, DataFrame)
        JuMP.@constraint(model,
        inv_gamma[a in 1:E],
        gamma[a] == inv[a, :Upgrade_Lvl]
        )
    end
    if storage_investments !== nothing
        inv = CSV.read(storage_investments, DataFrame)
        JuMP.@constraint(model,
        inv_s_power[i in 1:N],
        s_power[i] == inv[i, :Storage_Power]
        )
        JuMP.@constraint(model,
        inv_s_energy[i in 1:N],
        s_energy[i] == inv[i, :Storage_Energy]
        )
        JuMP.@constraint(model,
        inv_sigma[i in 1:N],
        sigma[i] == ((inv[i, :Storage_Power] != 0 || inv[i, :Storage_Energy] != 0) ? 1 : 0)
        )
    end

    # global power balance
    JuMP.@constraint(model, 
        power_balance[r in 1:R, t in 1:T],
        sum(pg[r,g,t] for g in 1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i in 1:N)
        + sum(ue[r,i,t] for i in 1:N) == 0
    )

    # flow constraints, F = PTDF * (Pg - Pd + U)
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    JuMP.@constraint(model, 
    flow_lb[r in 1:R, a in rate_a_nonzero, t in 1:T],
    sum(ptdf[a, i - 1] * (
        sum(pg[r, g, t] for g in [gen_id for gen_type in values(data["bus"]["$i"]["gen"]) for gen_id in gen_type]) 
        - data["bus"]["$i"]["load"]["$r"][t] 
        + ue[r, i, t]) for i in 2:N)
        >= -data["branch"]["$a"]["rate_a"] - gamma[a] * get_capacity_increment(data, a)
    )

    JuMP.@constraint(model, 
        flow_ub[r in 1:R, a in rate_a_nonzero, t in 1:T],
        sum(ptdf[a, i - 1] * (
            sum(pg[r, g, t] for g in [gen_id for gen_type in values(data["bus"]["$i"]["gen"]) for gen_id in gen_type]) 
            - data["bus"]["$i"]["load"]["$r"][t] 
            + ue[r, i, t]) for i in 2:N)
        <= data["branch"]["$a"]["rate_a"] + gamma[a] * get_capacity_increment(data, a)
    )

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    #
    #   III. Objective
    #

    operational_weight = 1
    if haskey(data["param"], "operational_weight")
        operational_weight = data["param"]["operational_weight"]
    end
    
    # Define the base objective expression
    base_objective = (sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
    sum(sigma[i] for i in 1:N) * data["param"]["storage_fixed_cost"] + 
    sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
    sum(
        data["param"]["representative_prob"][r] *
        (
            sum(
                sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g in 1:G) +
                sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i in 1:N)
            for t in 1:T)
        )
    for r in 1:R) * operational_weight)

    # Set the objective once
    JuMP.@objective(model, Min, base_objective)

    return model
end