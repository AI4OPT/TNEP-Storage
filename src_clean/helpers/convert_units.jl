# Function to convert cost coefficients to per unit
function convert_to_pu(cost, ncost, base_power = 100)
    cost_pu = []

    for i in 1:ncost
        if i == 1
            push!(cost_pu, cost[i])  # a term (constant)
        else
            push!(cost_pu, cost[i] * base_power^(i-1))  # scale by base power
        end
    end
    
    return cost_pu
end

"""
This function converts the powermodels file into radians and p.u.
"""
function convert_units(data, base_power  = 100)
    MW_TO_PU = 1 / base_power
    PU_TO_MW = base_power
    # first convert the load profiles into p.u.
    for (bus_key, bus_value) in data["bus"]
        for (load_key, load_value) in bus_value["load"]
            data["bus"][bus_key]["load"][load_key] = load_value .* MW_TO_PU
        end
    end

    # convert rate_a of each branch into p.u.
    for (branch_key, branch_value) in data["branch"]
        branch_value["rate_a"] *= MW_TO_PU
    end

    renewable_types = data["param"]["renewable_types"]
    nonrenewable_types = data["param"]["nonrenewable_types"]
    # convert all pmax and pmin into p.u.
    # convert renewable profiles or foreign imports into p.u.
    # convert cost curves of nonrenewable profiles into p.u.
    for (gen_key, gen_value) in data["gen"]
        gen_value["pmax"] *= MW_TO_PU
        gen_value["pmin"] *= MW_TO_PU
        if gen_value["gen_type"] in renewable_types || gen_value["gen_type"] == "foreign"
            for (profile_key, profile_value) in gen_value["profile"]
                data["gen"][gen_key]["profile"][profile_key] = profile_value .* MW_TO_PU
            end
        else
            data["gen"][gen_key]["cost"] = convert_to_pu(gen_value["cost"], gen_value["ncost"])
        end
    end

    return data
end


