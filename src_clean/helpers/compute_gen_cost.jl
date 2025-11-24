function compute_gen_cost(pg, gen_data)
    model = gen_data["model"]
    ncost = gen_data["ncost"]
    cost = gen_data["cost"]
    pmax = gen_data["pmax"]
    pmin = gen_data["pmin"]
    
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

        # return sum([cost[k] * pg^(k-1) for k in 1:ncost])

        """
        # this will return an over-approximation of the linearization of the quadratic cost curve
        if pmax == pmin
            return 0
        end
        a = cost[1]
        b = cost[2]
        c = cost[3]
        C_Pmin = a + b * pmin + c * pmin^2
        C_Pmax = a + b * pmax + c * pmax^2
        m = (C_Pmax - C_Pmin) / (pmax - pmin)
       return C_Pmin + m * (pg - pmin)"""

       # this is the correct linear version (deleting the quadratic term)
       # return sum([cost[k] * pg^(k-1) for k in 1:2])
       # this is the linear version that ignores the constant cost in the objective
       return cost[2] * pg
    else
        error("Unsupported cost model type: $model")
    end
end


        