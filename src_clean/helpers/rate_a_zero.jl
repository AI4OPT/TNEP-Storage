function get_rate_a_zero(data)
    rate_a_zero = Set{String}()
    rate_a_nonzero = Set{String}()
    for (key, value) in data["branch"]
        if value["rate_a"] == 0
            push!(rate_a_zero, key)
        else
            push!(rate_a_nonzero, key)
        end
    end

    return rate_a_zero, rate_a_nonzero
end

function get_capacity_increment(data, arc)
    cap_upgrade_increment = data["param"]["cap_upgrade_increment"]
    if haskey(data["param"], "cap_percent") && data["param"]["cap_percent"] == true
        cap_upgrade_increment = data["param"]["cap_upgrade_increment"] * data["branch"]["$(arc)"]["rate_a"]
    end

    return cap_upgrade_increment
end

