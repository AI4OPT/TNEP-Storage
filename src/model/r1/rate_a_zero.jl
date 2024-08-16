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

