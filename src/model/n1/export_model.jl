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
    num_nodes = length(data["bus"])
    num_hours = data["param"]["num_hours"]

    # Retrieve model values once to avoid multiple calls to value()
    ue = value.(model[:ue])[rep_index, :, :]
    ch = value.(model[:ch])[rep_index, :, :]

    # Check if 'oe' exists in the model
    has_oe = haskey(model, :oe)
    oe = has_oe ? value.(model[:oe])[rep_index, :, :] : zeros(size(ue))

    # Check if 'dis' exists in the model
    has_dis = haskey(model, :dis)
    dis = has_dis ? value.(model[:dis])[rep_index, :, :] : zeros(size(ch))

    # Infer 'dis' from negative 'ch' values if it does not exist
    if !has_dis
        dis .= max.(-ch, 0)  # Set dis to -ch when ch < 0, else 0
        ch .= max.(ch, 0)     # Keep only positive charge values
    end

    ### **Optimized Power Flow Retrieval**
    has_pf = haskey(model, :pf)
    if has_pf
        # Use precomputed power flows (DCOPF case)
        pfs = value.(model[:pf])[rep_index, :, :]
        # Since nodal balance is already enforced in DCOPF, we can directly use ue - oe
        imbalance = ue - oe
    else
        # Compute power flows using PTDF
        pg_values = value.(model[:pg])[rep_index:rep_index, :, :]
        ue_values = value.(model[:ue])[rep_index:rep_index, :, :]
        ch_values = value.(model[:ch])[rep_index:rep_index, :, :]
        
        num_branches = length(data["branch"])
        flows = zeros(num_branches, 1, num_hours)
        
        compute_flows!(flows, pg_values, ue_values, ch_values, data, model.ext[:PTDF])
        pfs = flows[:, 1, :]

        ### **Precompute Generator Dispatch**
        pg_values = value.(model[:pg])[rep_index, :, :]

        ### **Precompute Nodal Connectivity for Efficient Lookup**
        arcs_from = Dict(i => [a for a in get(data["arcs_from"], "$i", []) if data["branch"]["$a"]["f_bus"] == i] for i in 1:num_nodes)
        arcs_to = Dict(i => [a for a in get(data["arcs_from"], "$i", []) if data["branch"]["$a"]["t_bus"] == i] for i in 1:num_nodes)

        ### **Compute Nodal Power Balance Efficiently (PTDF Case Only)**
        imbalance = zeros(num_nodes, num_hours)
        for bus in 1:num_nodes, hour in 1:num_hours
            generation = sum(pg_values[g, hour] for g in 1:length(data["gen"]) if data["gen"]["$g"]["gen_bus"] == bus; init=0.0)
            demand = data["bus"]["$bus"]["load"]["$rep_index"][hour]
            discharge = dis[bus, hour]
            charge = ch[bus, hour]

            flow_in = sum(pfs[branch, hour] for branch in get(arcs_to, bus, []); init=0.0)
            flow_out = sum(pfs[branch, hour] for branch in get(arcs_from, bus, []); init=0.0)

            imbalance[bus, hour] = (generation + discharge + ue[bus, hour] + flow_in) - (flow_out + demand + oe[bus, hour] + charge)
        end
    end

    ### **Data Mapping for Export**
    node_indices_mapping = ["$i" for i in 1:num_nodes]
    node_names_mapping = [data["bus"]["$i"]["bus_name"] for i in 1:num_nodes]
    lat_mapping = [data["bus"]["$i"]["lat"] for i in 1:num_nodes]
    lon_mapping = [data["bus"]["$i"]["lon"] for i in 1:num_nodes]
    hours_mapping = 1:num_hours

    ### **Column Names**
    column_names = if has_oe
        [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance, :Over_Energy, :Under_Energy, :Charge, :Discharge]
    else
        [:Node_Index, :Node_Name, :Lat, :Lon, :Hour, :Energy_Imbalance, :Charge, :Discharge]
    end

    ### **Flatten Data Efficiently**
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

    ### **Generator Output by Fuel Type**
    pg_values = value.(model[:pg])[rep_index, :, :]
    gen_types = vcat(data["param"]["renewable_types"], data["param"]["nonrenewable_types"])
    gen_array = [vec(transpose([gen_energy(bus, gen_type, hour, pg_values, data) for bus in 1:num_nodes, hour in 1:num_hours])) for gen_type in gen_types]
    df_generation = DataFrame(hcat(gen_array...), Symbol.(gen_types))

    ### **Renewable Production by Bus**
    dfs = [get_renewable_production_by_bus("$i", data, rep_index) for i in 1:num_nodes]
    df_renewable_prod = vcat(dfs...)
    rename_columns_with_suffix!(df_renewable_prod, "_production")

    ### **Combine & Save CSV**
    df_energy = hcat(df_imbalance, df_generation, df_renewable_prod)
    df_energy = round_df(df_energy)
    
    datestring = data["param"]["dates"][rep_index]
    CSV.write(joinpath(simdir, "output", datestring, "energy.csv"), df_energy)
end

function export_investments_csv(data, gammas, s_energys;
    s_powers=nothing,
    output_dir="output",
    file_suffix="")
    
    # Create the output directory if it doesn't exist
    mkpath(output_dir)
    
    # Construct filenames with optional suffix
    lines_filename = isempty(file_suffix) ? "line_investments.csv" : "line_investments_$(file_suffix).csv"
    storage_filename = isempty(file_suffix) ? "storage_investments.csv" : "storage_investments_$(file_suffix).csv"
    
    df_lines = DataFrame(
        Branch_Index = ["$i" for i in 1:length(data["branch"])],
        Lat1 = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lat"] for i in 1:length(data["branch"])],
        Lon1 = [data["bus"][string(data["branch"]["$i"]["f_bus"])]["lon"] for i in 1:length(data["branch"])],
        Lat2 = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lat"] for i in 1:length(data["branch"])],
        Lon2 = [data["bus"][string(data["branch"]["$i"]["t_bus"])]["lon"] for i in 1:length(data["branch"])],
        Rate_A = [data["branch"]["$i"]["rate_a"] for i in 1:length(data["branch"])],
        Upgrade_Lvl = [gammas[i] for i in 1:length(data["branch"])]
    )
    
    # DON'T ROUND THE LINES ANYMORE
    # df_lines = round_df(df_lines)
    CSV.write(joinpath(output_dir, lines_filename), df_lines)
    
    # Create storage DataFrame - only include Storage_Power column if s_powers is provided
    df_storage_dict = Dict(
        :Node_Index => ["$i" for i in 1:length(data["bus"])],
        :Node_Name => [data["bus"]["$i"]["bus_name"] for i in 1:length(data["bus"])],
        :Lat => [data["bus"]["$i"]["lat"] for i in 1:length(data["bus"])],
        :Lon => [data["bus"]["$i"]["lon"] for i in 1:length(data["bus"])],
        :Storage_Energy => [s_energys[i] for i in 1:length(data["bus"])]
    )
    
    # Add Storage_Power column only if s_powers is provided
    if s_powers !== nothing
        df_storage_dict[:Storage_Power] = [s_powers[i] for i in 1:length(data["bus"])]
    end
    
    df_storage = DataFrame(df_storage_dict)
    
    # DON'T ROUND THE STORAGE
    # df_storage = round_df(df_storage)
    CSV.write(joinpath(output_dir, storage_filename), df_storage)
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

    # DON'T ROUND FLOW ANYMORE
    # df_flow = round_df(df_flow)
    datestring = data["param"]["dates"][rep_index]
    CSV.write(joinpath(simdir, "output", datestring, "flow.csv"), df_flow)
end

function export_model(simdir, model, data)
    # Check if s_power exists in the model and call export_investments_csv accordingly
    if haskey(model, :s_power)
        export_investments_csv(data, value.(model[:gamma]), value.(model[:s_energy]), 
                              s_powers=value.(model[:s_power]), 
                              output_dir=joinpath(simdir, "output"))
    else
        export_investments_csv(data, value.(model[:gamma]), value.(model[:s_energy]), 
                              output_dir=joinpath(simdir, "output"))
    end
    
    for rep_index in 1:length(data["param"]["dates"])
        export_energy_csv(simdir, model, data, rep_index)
        export_flow(simdir, model, data, rep_index)
    end
end
