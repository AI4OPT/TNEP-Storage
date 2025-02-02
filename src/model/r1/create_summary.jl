include("storage_candidates/naive_candidates.jl")

function write_summary_to_csv(simdir, model, data)

    # Get number of line investments
    df_line_investments = CSV.read(joinpath(simdir, "output", "line_investments.csv"), DataFrame)
    num_line_investments = count(!=(0), df_line_investments[!, :Upgrade_Lvl])

    # Get number of storage investments
    df_storage_investments = CSV.read(joinpath(simdir, "output", "storage_investments.csv"), DataFrame)
    num_storage_investments = count(!=(0), df_storage_investments[!, :Storage_Energy])

    # Total storage capacity
    total_storage_capacity = sum(df_storage_investments[!, :Storage_Energy])

    # Average storage capacity
    avg_storage_capacity = total_storage_capacity / num_storage_investments

    # Get line investment costs
    col_data = df_line_investments[!, :Upgrade_Lvl]
    E = length(col_data)
    line_investment_costs = sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * col_data[a] for a in 1:E)

    # Get storage investment costs
    power_data = df_storage_investments[!, :Storage_Power]
    energy_data = df_storage_investments[!, :Storage_Energy]
    N = length(energy_data)
    storage_investment_costs = (
        sum(power_data[i] * data["param"]["bess_power_cost"] 
        + energy_data[i] * data["param"]["bess_energy_cost"] for i in 1:N) 
        + num_storage_investments * data["param"]["storage_fixed_cost"]
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
    
    # Calculate under-served penalty
    under_served_penalty = operational_weight * sum(data["param"]["representative_prob"][r] * 
        sum(
            sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
        for t in 1:T)
    for r in 1:R)

    # Number of storage candidates, if available
    storage_cand_count = length(data["bus"])
    if haskey(data["param"], "candidate_no_upgrades_dir")
        bus_set = intersect_storage_candidates(data, data["param"]["candidate_no_upgrades_dir"])
        storage_cand_count = length(bus_set)
    end

    # Create base summary data
    summary_data_dict = Dict(
        "under_served_penalty" => under_served_penalty,
        "generation_costs" => generation_costs,
        "storage_investment_costs" => storage_investment_costs,
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
            "storage_investment_costs", "line_investment_costs", "num_line_investments",
            "num_storage_investments", "total_storage_capacity", "avg_storage_capacity",
            "storage_cand_count"]
            
    for var in order
        if haskey(summary_data_dict, var)
            push!(summary_data, (var, summary_data_dict[var]))
        end
    end

    CSV.write(joinpath(simdir, "output", "summary_data.csv"), summary_data)
end