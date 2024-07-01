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

function get_renewable_production_by_bus(bus, data)
    renewable_production = Dict()
    renewable_types = data["param"]["renewable_types"]
    bus_gens = data["bus"][bus]["gen"]

    for renewable_type in renewable_types
        renewable_production[renewable_type] = zeros(data["param"]["num_hours"])
        if haskey(bus_gens, renewable_type)
            for gen_id in bus_gens[renewable_type]
                profile = data["gen"]["$gen_id"]["profile"]
                avg_time_series = mean(reduce(hcat, values(profile)), dims=2)
                renewable_production[renewable_type] += avg_time_series
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

function export_energy_csv(simdir, model, data)

    ue = mean(value.(model[:ue]), dims=1)
    oe = mean(value.(model[:oe]), dims=1)
    ch = mean(value.(model[:ch]), dims=1)
    dis = mean(value.(model[:dis]), dims=1)
    # ch = mean(value.(model[:ch]), dims=1) * data["param"]["bess_efficiency"]
    # dis = mean(value.(model[:dis]), dims=1) / data["param"]["bess_efficiency"]
    imbalance = oe - ue # 1 x num_nodes x num_hours

    num_nodes = size(imbalance, 2)
    num_hours = size(imbalance, 3)

    # Mappings (you need to define these based on your data)
    node_indices_mapping = ["$i" for i in 1:length(data["bus"])]
    node_names_mapping = [data["bus"]["$i"]["bus_name"] for i in 1:length(data["bus"])]
    lat_mapping = [data["bus"]["$i"]["lat"] for i in 1:length(data["bus"])]
    lon_mapping = [data["bus"]["$i"]["lon"] for i in 1:length(data["bus"])]
    hours_mapping = 1:data["param"]["num_hours"]

    column_names = [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance, :Charge, :Discharge]

    # Flatten the array and create DataFrame
    energy = [(node_indices_mapping[bus], 
        node_names_mapping[bus],
        lat_mapping[bus],
        lon_mapping[bus],
        hours_mapping[hour], 
        imbalance[1,bus,hour],
        ch[1,bus,hour],
        dis[1,bus,hour]) 
            for bus in 1:num_nodes, hour in 1:num_hours]

    energy_flattened = [vec for row in eachrow(energy) for vec in row]
    df_imbalance = DataFrame(energy_flattened, Symbol.(column_names))


    # Now create the generator output information by fuel type
    pg = mean(value.(model[:pg]), dims=1) # 1 x 210 x 24
    
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

    # Now create the renewable production information by fuel type (before curtailment)
    dfs = DataFrame[]
    for i in 1:length(data["bus"])
        df_bus = get_renewable_production_by_bus("$i", data)
        push!(dfs, df_bus)
    end
    df_renewable_prod = vcat(dfs...)
    rename_columns_with_suffix!(df_renewable_prod, "_production")

    # Horizontally concatenate the imbalance and generation information
    df_energy = hcat(df_imbalance, df_generation, df_renewable_prod) 
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

function export_flow(simdir, model, data)
    pfs = mean(value.(model[:pf]), dims=1)
    num_branches = size(pfs, 2)
    num_hours = size(pfs, 3)

    # Mappings (you need to define these based on your data)
    branch_indices_mapping = ["$i" for i in 1:length(data["branch"])]
    lat1_mapping = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lat"] for i in 1:length(data["branch"])]
    lon1_mapping = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lon"] for i in 1:length(data["branch"])]
    lat2_mapping = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lat"] for i in 1:length(data["branch"])]
    lon2_mapping = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lon"] for i in 1:length(data["branch"])]
    thermal_limits_mapping = [data["branch"]["$i"]["rate_a"] for i in 1:length(data["branch"])]

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
        pfs[1,branch,hour]) 
            for branch in 1:num_branches, hour in 1:num_hours]

    flow_flattened = [vec for row in eachrow(flow) for vec in row]
    df_flow = DataFrame(flow_flattened, Symbol.(column_names))

    df_flow = round_df(df_flow)
    CSV.write(joinpath(simdir, "output", "flow.csv"), df_flow)
end

# ue, oe, pf, dis, ch = value.(model[:ue]), value.(model[:oe]), value.(model[:pf]), value.(model[:dis]), value.(model[:ch])
