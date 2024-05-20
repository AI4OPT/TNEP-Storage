function compute_gen_cost(pg, gen_data)
    model = gen_data["model"]
    ncost = gen_data["ncost"]
    cost = gen_data["cost"]
    
    if model == 1
        # Piecewise Linear Model
        # Assuming cost data is sorted by power output breakpoints
        # and pg is within the bounds of the cost data.
        for i in 2:ncost
            if pg <= cost[i][1]
                # Linear interpolation between two points
                m = (cost[i][2] - cost[i-1][2]) / (cost[i][1] - cost[i-1][1])
                return cost[i-1][2] + m * (pg - cost[i-1][1])
            end
        end
    elseif model == 2
        # Polynomial Cost Model
        # The cost array contains coefficients for the polynomial
        # where the highest degree coefficient comes first.
        return sum([cost[k] * pg^(k-1) for k in 1:ncost])
    else
        error("Unsupported cost model type: $model")
    end
end
        