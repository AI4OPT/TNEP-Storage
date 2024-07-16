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
    oe = value.(model[:oe])
    R = data["param"]["num_representatives"]
    G = length(data["gen"])
    T = data["param"]["num_hours"]
    N = length(data["bus"])

    # Generation costs
    generation_costs = sum(data["representative_prob"][r] * 
        sum(
            sum(compute_gen_cost(pg[r,g,t], data["gen"]["$g"]) for g in 1:G)
        for t in 1:T)
    for r in 1:R)

    # Overgeneration penalty
    over_generation_penalty = sum(data["representative_prob"][r] * 
        sum(
            sum(data["param"]["over_generated_penalty"] * oe[r,i,t] for i in 1:N)
        for t in 1:T)
    for r in 1:R)
    
    # Underserved penalty
    under_served_penalty = sum(data["representative_prob"][r] * 
        sum(
            sum(data["param"]["under_served_penalty"] * ue[r,i,t] for i in 1:N)
        for t in 1:T)
    for r in 1:R)

    summary_data = DataFrame(
        Variable = [
            "over_generation_penalty", 
            "under_served_penalty", 
            "generation_costs", 
            "storage_investment_costs", 
            "line_investment_costs", 
            "num_line_investments", 
            "num_storage_investments", 
            "total_storage_capacity", 
            "avg_storage_capacity"
        ],
        Value = [
            over_generation_penalty, 
            under_served_penalty, 
            generation_costs, 
            storage_investment_costs, 
            line_investment_costs, 
            num_line_investments, 
            num_storage_investments, 
            total_storage_capacity, 
            avg_storage_capacity
        ]
    )

    CSV.write(joinpath(simdir, "output", "summary_data.csv"), summary_data)
end
