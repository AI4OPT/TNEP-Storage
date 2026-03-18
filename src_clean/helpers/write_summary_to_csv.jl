function get_capacity_increment(data, arc)
    cap_upgrade_increment = data["param"]["cap_upgrade_increment"]
    if haskey(data["param"], "cap_percent") && data["param"]["cap_percent"] == true
        cap_upgrade_increment = data["param"]["cap_upgrade_increment"] * data["branch"]["$(arc)"]["rate_a"]
    end

    return cap_upgrade_increment
end

function write_summary_to_csv(simdir, model, data)

    # Get number of line investments
    df_line_investments = CSV.read(joinpath(simdir, "output", "line_investments.csv"), DataFrame)
    num_line_investments = count(x -> abs(x) > 1e-5, df_line_investments[!, :Upgrade_Lvl])

    # Get number of storage investments
    df_storage_investments = CSV.read(joinpath(simdir, "output", "storage_investments.csv"), DataFrame)
    num_storage_investments = count(x -> abs(x) > 1e-5, df_storage_investments[!, :Storage_Energy])

    # Total storage capacity
    total_storage_capacity = sum(df_storage_investments[!, :Storage_Energy])

    # Average storage capacity
    avg_storage_capacity = total_storage_capacity / num_storage_investments

    # Get line investment costs
    col_data = df_line_investments[!, :Upgrade_Lvl]
    E = length(col_data)
    line_investment_costs = sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * data["branch"]["$a"]["distance"] * col_data[a] for a in 1:E)

    # Get storage investment costs
    energy_data = df_storage_investments[!, :Storage_Energy]
    N = length(energy_data)
    storage_investment_costs = (
        sum(energy_data[i] * data["param"]["bess_energy_cost"] for i in 1:N) * data["param"]["storage_energy_size"] 
    )

    # Model variable values
    pg = value.(model[:pg])
    ue = value.(model[:ue])
    R = data["param"]["num_representatives"]
    G = length(data["gen"])
    T = data["param"]["num_hours"]
    N = length(data["bus"])

    operational_weight = get(data["param"], "operational_weight", 1)

    # Generation costs
    generation_costs = operational_weight * sum(data["param"]["representative_prob"][r] * 
        sum(
            sum(compute_gen_cost(pg[r,g,t], data["gen"]["$g"]) for g in 1:G)
        for t in 1:T)
    for r in 1:R)

    # Check if over-energy exists in the model
    has_oe = haskey(model, :oe)
    
    # Initialize penalties
    over_generation_penalty = 0.0
    under_served_penalty = 0.0
    
    if has_oe
        # Calculate over-generation penalty if oe exists
        oe = value.(model[:oe])
        over_generation_penalty = operational_weight * sum(data["param"]["representative_prob"][r] * 
            sum(
                sum(data["param"]["over_generated_penalty"] * oe[r,i,t] for i in 1:N)
            for t in 1:T)
        for r in 1:R)
    end


    ch = value.(model[:ch])
    dis = value.(model[:dis])
    # Storage operation costs
    storage_operation_costs = operational_weight * sum(data["param"]["representative_prob"][r] * 
            sum(
                sum(get(data["param"], "storage_operation_cost", 0.0) * (ch[r,i,t] + dis[r,i,t]) for i=1:N)
            for t in 1:T)
        for r in 1:R)

    
    # Calculate under-served penalty
    under_served_penalty = operational_weight * sum(data["param"]["representative_prob"][r] * 
        sum(
            sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
        for t in 1:T)
    for r in 1:R)

    # Number of storage candidates, if available
    storage_cand_count = length(data["bus"])

    # Create base summary data
    summary_data_dict = Dict(
        "under_served_penalty" => under_served_penalty,
        "generation_costs" => generation_costs,
        "storage_investment_costs" => storage_investment_costs,
        "storage_operation_costs" => storage_operation_costs,
        "line_investment_costs" => line_investment_costs,
        "num_line_investments" => num_line_investments,
        "num_storage_investments" => num_storage_investments,
        "total_storage_capacity" => total_storage_capacity,
        "avg_storage_capacity" => avg_storage_capacity,
        "storage_cand_count" => storage_cand_count
    )
    
    # Add over_generation_penalty only if oe exists
    if has_oe
        summary_data_dict["over_generation_penalty"] = over_generation_penalty
    end

    # Convert to DataFrame
    summary_data = DataFrame(
        Variable = String[],
        Value = Float64[]
    )
    
    # Add rows in desired order
    order = ["over_generation_penalty", "under_served_penalty", "generation_costs", 
            "storage_investment_costs", "storage_operation_costs", "line_investment_costs", "num_line_investments",
            "num_storage_investments", "total_storage_capacity", "avg_storage_capacity",
            "storage_cand_count"]
            
    for var in order
        if haskey(summary_data_dict, var)
            push!(summary_data, (var, summary_data_dict[var]))
        end
    end

    CSV.write(joinpath(simdir, "output", "summary_data.csv"), summary_data)
end