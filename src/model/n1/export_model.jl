using CSV
using Statistics
using DataFrames

function round_df(df)
    for col_name in names(df)
        # Check if the column type is a subtype of AbstractFloat and the column name is not Lat or Lon
        if eltype(df[!, col_name]) <: AbstractFloat
            df[!, col_name] = round.(df[!, col_name], digits=3)
        end
    end

    return df
end

function gen_energy(bus, gen_type, hour, pg, data)
    if !haskey(data["bus"]["$bus"]["gen"], gen_type)
        return 0.0
    else
        return sum(pg[gen_id,hour] for gen_id in data["bus"]["$bus"]["gen"][gen_type])
    end
end

function get_renewable_production_by_bus(bus, data, rep_index)
    renewable_production = Dict()
    renewable_types = data["param"]["renewable_types"]
    bus_gens = data["bus"][bus]["gen"]

    for renewable_type in renewable_types
        renewable_production[renewable_type] = zeros(data["param"]["num_hours"])
        if haskey(bus_gens, renewable_type)
            for gen_id in bus_gens[renewable_type]
                profile = data["gen"]["$gen_id"]["profile"]["$rep_index"]
                renewable_production[renewable_type] += profile
            end
        end
    end
    renewable_production = Dict(k => vec(v) for (k, v) in renewable_production)
    df = DataFrame(renewable_production)
    df = df[:, renewable_types]
    return df
end

function rename_columns_with_suffix!(df::DataFrame, suffix::String)
    old_names = names(df)
    new_names = Symbol[]
    for name in old_names
        push!(new_names, Symbol(string(name) * suffix))
    end
    rename!(df, new_names)
end

function export_energy_csv(simdir, model, data, rep_index)
    ue = value.(model[:ue])[rep_index, :, :]
    ch = value.(model[:ch])[rep_index, :, :]
    
    # Check if 'oe' exists in the model
    has_oe = haskey(model, :oe)
    oe = has_oe ? value.(model[:oe])[rep_index, :, :] : zeros(size(ue))
    
    # Check if 'dis' exists in the model
    has_dis = haskey(model, :dis)
    dis = has_dis ? value.(model[:dis])[rep_index, :, :] : zeros(size(ch))

    # If 'dis' doesn't exist, infer discharge from negative 'ch' values
    if !has_dis
        for bus in 1:size(ch, 1), hour in 1:size(ch, 2)
            if ch[bus, hour] < 0
                dis[bus, hour] = -ch[bus, hour]  # Assign negative charge to discharge
                ch[bus, hour] = 0                # Set charge to zero
            end
        end
    end

    imbalance = oe - ue  # num_nodes x num_hours

    num_nodes = size(imbalance, 1)
    num_hours = size(imbalance, 2)

    # Mappings
    node_indices_mapping = ["$i" for i in 1:length(data["bus"])]
    node_names_mapping = [data["bus"]["$i"]["bus_name"] for i in 1:length(data["bus"])]
    lat_mapping = [data["bus"]["$i"]["lat"] for i in 1:length(data["bus"])]
    lon_mapping = [data["bus"]["$i"]["lon"] for i in 1:length(data["bus"])]
    hours_mapping = 1:data["param"]["num_hours"]

    # Add OE to column names if it exists
    column_names = if has_oe
        [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance, :Over_Energy, :Under_Energy, :Charge, :Discharge]
    else
        [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance, :Charge, :Discharge]
    end

    # Flatten the array and create DataFrame with optional OE
    energy = if has_oe
        [(node_indices_mapping[bus],
          node_names_mapping[bus],
          lat_mapping[bus],
          lon_mapping[bus],
          hours_mapping[hour],
          imbalance[bus, hour],
          oe[bus, hour],
          ue[bus, hour],
          ch[bus, hour],
          dis[bus, hour])
          for bus in 1:num_nodes, hour in 1:num_hours]
    else
        [(node_indices_mapping[bus],
          node_names_mapping[bus],
          lat_mapping[bus],
          lon_mapping[bus],
          hours_mapping[hour],
          imbalance[bus, hour],
          ch[bus, hour],
          dis[bus, hour])
          for bus in 1:num_nodes, hour in 1:num_hours]
    end

    energy_flattened = [vec for row in eachrow(energy) for vec in row]
    df_imbalance = DataFrame(energy_flattened, Symbol.(column_names))

    # Rest of the function remains the same
    # Generator output by fuel type
    pg = value.(model[:pg])[rep_index, :, :]
    gen_types = vcat(data["param"]["renewable_types"], data["param"]["nonrenewable_types"])
    gen_array = Float64[]

    for gen_type in gen_types
        x = vec(transpose([gen_energy(bus, gen_type, hour, pg, data) for bus in 1:num_nodes, hour in 1:num_hours]))
        if isempty(gen_array)
            gen_array = x
        else
            gen_array = hcat(gen_array, x)
        end
    end
    df_generation = DataFrame(gen_array, Symbol.(gen_types))

    # Renewable production by bus
    dfs = DataFrame[]
    for i in 1:length(data["bus"])
        df_bus = get_renewable_production_by_bus("$i", data, rep_index)
        push!(dfs, df_bus)
    end
    df_renewable_prod = vcat(dfs...)
    rename_columns_with_suffix!(df_renewable_prod, "_production")

    # Combine all DataFrames
    df_energy = hcat(df_imbalance, df_generation, df_renewable_prod)
    df_energy = round_df(df_energy)
    
    datestring = data["param"]["dates"][rep_index]
    CSV.write(joinpath(simdir, "output", datestring, "energy.csv"), df_energy)
end

function export_investments_csv(simdir, model, data)
    gammas = value.(model[:gamma])
    s_powers = value.(model[:s_power])
    s_energys = value.(model[:s_energy])

    df_lines = DataFrame(
        Branch_Index = ["$i" for i in 1:length(data["branch"])],
        Lat1 = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lat"] for i in 1:length(data["branch"])],
        Lon1 = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lon"] for i in 1:length(data["branch"])],
        Lat2 = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lat"] for i in 1:length(data["branch"])],
        Lon2 = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lon"] for i in 1:length(data["branch"])],
        Rate_A = [data["branch"]["$i"]["rate_a"] for i in 1:length(data["branch"])],
        Upgrade_Lvl = [gammas[i] for i in 1:length(data["branch"])]
    )

    df_lines = round_df(df_lines)
    CSV.write(joinpath(simdir, "output", "line_investments.csv"), df_lines)

    df_storage = DataFrame(
        Node_Index = ["$i" for i in 1:length(data["bus"])],
        Node_Name = [data["bus"]["$i"]["bus_name"] for i in 1:length(data["bus"])],
        Lat = [data["bus"]["$i"]["lat"] for i in 1:length(data["bus"])],
        Lon = [data["bus"]["$i"]["lon"] for i in 1:length(data["bus"])],
        Storage_Power = [s_powers[i] for i in 1:length(data["bus"])],
        Storage_Energy = [s_energys[i] for i in 1:length(data["bus"])]
    )

    df_storage = round_df(df_storage)
    CSV.write(joinpath(simdir, "output", "storage_investments.csv"), df_storage)
end

function export_flow(simdir, model, data, rep_index)
    # Check if we're using PTDF model or direct flow variables
    if haskey(model, :pf)
        # Original model with direct flow variables
        pfs = value.(model[:pf])[rep_index, :, :]
    else
        # PTDF model - calculate flows
        pg_values = value.(model[:pg])[rep_index:rep_index, :, :]
        ue_values = value.(model[:ue])[rep_index:rep_index, :, :]
        ch_values = value.(model[:ch])[rep_index:rep_index, :, :]
        
        num_branches = length(data["branch"])
        num_hours = data["param"]["num_hours"]
        flows = zeros(num_branches, 1, num_hours)
        
        compute_flows!(flows, pg_values, ue_values, ch_values, data, model.ext[:PTDF])
        pfs = flows[:, 1, :]
    end

    num_branches = size(pfs, 1)
    num_hours = size(pfs, 2)

    # Mappings (same as before)
    branch_indices_mapping = ["$i" for i in 1:length(data["branch"])]
    lat1_mapping = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lat"] for i in 1:length(data["branch"])]
    lon1_mapping = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lon"] for i in 1:length(data["branch"])]
    lat2_mapping = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lat"] for i in 1:length(data["branch"])]
    lon2_mapping = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lon"] for i in 1:length(data["branch"])]
    
    # Get thermal limits including any upgrades if they exist
    if haskey(model, :gamma)
        gamma_values = value.(model[:gamma])
        thermal_limits_mapping = [
            data["branch"]["$i"]["rate_a"] + gamma_values[i] * get_capacity_increment(data, i) 
            for i in 1:length(data["branch"])
        ]
    else
        thermal_limits_mapping = [data["branch"]["$i"]["rate_a"] for i in 1:length(data["branch"])]
    end

    hours_mapping = 1:data["param"]["num_hours"]

    column_names = [:Branch_Index, :Lat1, :Lon1, :Lat2, :Lon2, :Hour, :Rate_A, :Power_Flow]

    # Flatten the array and create DataFrame
    flow = [(branch_indices_mapping[branch], 
        lat1_mapping[branch],
        lon1_mapping[branch],
        lat2_mapping[branch],
        lon2_mapping[branch],
        hours_mapping[hour], 
        thermal_limits_mapping[branch],
        pfs[branch,hour]) 
            for branch in 1:num_branches, hour in 1:num_hours]

    flow_flattened = [vec for row in eachrow(flow) for vec in row]
    df_flow = DataFrame(flow_flattened, Symbol.(column_names))

    df_flow = round_df(df_flow)
    datestring = data["param"]["dates"][rep_index]
    CSV.write(joinpath(simdir, "output", datestring, "flow.csv"), df_flow)
end

function export_model(simdir, model, data)
    export_investments_csv(simdir, model, data)

    for rep_index in 1:length(data["param"]["dates"])
        export_energy_csv(simdir, model, data, rep_index)
        export_flow(simdir, model, data, rep_index)
    end
end
