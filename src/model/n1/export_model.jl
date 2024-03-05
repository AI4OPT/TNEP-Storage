using CSV
using Statistics
using DataFrames

function round_df(df)
    for col_name in names(df)
        # Check if the column type is a subtype of AbstractFloat
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
        return sum(pg[1,gen_id,hour] for gen_id in data["bus"]["$bus"]["gen"][gen_type])
    end
end

function export_energy_csv(simdir, model, data)

    ue = mean(value.(model[:ue]), dims=1)
    oe = mean(value.(model[:oe]), dims=1)
    imbalance = oe - ue # 1 x num_nodes x num_hours

    num_nodes = size(imbalance, 2)
    num_hours = size(imbalance, 3)

    # Mappings (you need to define these based on your data)
    node_indices_mapping = ["$i" for i in 1:length(data["bus"])]
    node_names_mapping = [data["bus"]["$i"]["bus_name"] for i in 1:length(data["bus"])]
    lat_mapping = [data["bus"]["$i"]["lat"] for i in 1:length(data["bus"])]
    lon_mapping = [data["bus"]["$i"]["lon"] for i in 1:length(data["bus"])]
    hours_mapping = 1:data["param"]["num_hours"]

    column_names = [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance]

    # Flatten the array and create DataFrame
    energy = [(node_indices_mapping[bus], 
        node_names_mapping[bus],
        lat_mapping[bus],
        lon_mapping[bus],
        hours_mapping[hour], 
        imbalance[1,bus,hour]) 
            for bus in 1:num_nodes, hour in 1:num_hours]

    energy_flattened = [vec for row in eachrow(energy) for vec in row]
    df_imbalance = DataFrame(energy_flattened, Symbol.(column_names))


    # Now create the generator information by fuel type
    pg = mean(value.(model[:pg]), dims=1) # 1 x 210 x 24
    
    gen_types = data["param"]["gen_types"]
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

    # Horizontally concatenate the imbalance and generation information
    df_energy = hcat(df_imbalance, df_generation) 
    df_energy = round_df(df_energy)
    CSV.write(joinpath(simdir, "output", "energy.csv"), df_energy)
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

# ue, oe, pf, dis, ch = value.(model[:ue]), value.(model[:oe]), value.(model[:pf]), value.(model[:dis]), value.(model[:ch])
