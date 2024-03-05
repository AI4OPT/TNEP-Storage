using JSON
using CSV
using DataFrames
include("../../helpers/haversine_distance.jl")
include("../../helpers/safe_parse_to_float.jl")

const n1_balance_file = "data/topology/n1/EIA930_BALANCE_2019_Jan_Jun.csv"
const n1_interchange_file = "data/topology/n1/EIA930_INTERCHANGE_2019_Jan_Jun.csv"
const n1_coords_json = "data/topology/n1/xycoords.json"
const production_costs = Dict(
    "units" => "dollar/MW",
    "coal" => 36.66,
    "natural_gas" => 36.66,
    "nuclear" => 23.73,
    "hydro" => 40.00,
    "solar" => 0.0,
    "wind" => 0.0
)

function add_generators_from_df(df, gen_type, p_data)
    num_gens = length(p_data["gen"])
    id = num_gens + 1

    for row in eachrow(df)
        bus = row[Symbol("Balancing Authority")]
        bus_id = p_data["bus_name_to_id"][bus]
        capacity = row[2]  # Assuming the max capacity is the third column

        if capacity > 0
            gen_id = id # Create a unique generator ID
            # Add the generator
            p_data["gen"][gen_id] = Dict(
                "gen_bus" => bus_id, 
                "pmax" => capacity,
                "pmin" => 0.0,
                "gen_type" => gen_type,
                "model" => 1,
                "ncost" => 2,
                "cost"  => [0, production_costs[gen_type]],
                "profile" => Dict()
            )

            # Add generator to the bus information too
            p_data["bus"][bus_id]["gen"][gen_type] = [id]
            id += 1
        end
    end
end

function create_p_data(interchange_file, balance_file)
    p_data = Dict(
        "per_unit" => false,
        "baseMVA" => 100.0, 
        "multinetwork" => false,
        "name" => "US_wide_transport_model",
        "description" => "Low granularity representation of US electric grid structure; includes all regional operators",
        "source_type" => "EIA.gov",
        "source_version" => "2019_Jan_Jun",
        "gen" => Dict()
    )

    csv_file_path = interchange_file
    df = CSV.read(csv_file_path, DataFrame)

    # Get unique bus names
    bus_names = unique(df[:, "Balancing Authority"])

    # Initialize the branches dictionary and interchange data
    branches = Dict{String,Any}()
    branch_interchanges = Dict{String,Vector{Float64}}()

    # Create unique pairs and add branches with interchange data
    for row in eachrow(df)
        from_bus = row["Balancing Authority"]
        to_bus = row["Directly Interconnected Balancing Authority"]
        # Set default value to 0.0
        interchange_value = 0.0 
        # Check if the value is missing and parse if it's present
        interchange_value = safe_parse_to_float(row["Interchange (MW)"])
        if !ismissing(to_bus) && !ismissing(interchange_value) && from_bus != to_bus # Check for valid data
            # Create a sorted tuple to ensure uniqueness (A->B is same as B->A)
            pair = [String(bus) for bus in sort([from_bus, to_bus])]
            pair_key = pair[1] * "_" * pair[2]
            # If a pair is not in the dict yet, add it
            if pair[1] ∉ bus_names
                push!(bus_names, pair[1])
            end
            if pair[2] ∉ bus_names
                push!(bus_names, pair[2])
            end
            # Add interchange values to the interchange data array for the pair
            if haskey(branch_interchanges, pair_key)
                push!(branch_interchanges[pair_key], interchange_value)
            else
                branch_interchanges[pair_key] = [interchange_value]
            end
            # Create branch if it doesn't exist
            if !haskey(branches, pair_key)
                branches[pair_key] = Dict{Any,Any}("f_bus" => pair[1], "t_bus" => pair[2], "branch_id" => length(branches) + 1)
            end
        end
    end

    # Determine the maximum interchange for each branch to set rate_a
    for (pair_key, interchanges) in branch_interchanges
        # Calculate the maximum after taking the absolute value of each interchange
        max_interchange = maximum(abs.(interchanges))
        branches[pair_key]["rate_a"] = max_interchange # Set the maximum absolute interchange as rate_a for the branch
    end

    # Initialize the buses dictionary
    buses = Dict{Any,Any}()
    bus_name_to_id = Dict{Any,Any}()
    xy_coords = JSON.parsefile(v1_coords_json)
    for i in 1:length(bus_names)
        bus_name = bus_names[i]
        buses[i] = Dict("bus_type" => 1, 
            "bus_name" => bus_name, 
            "index" => i, 
            "lat" => xy_coords[bus_name][2], 
            "lon" => xy_coords[bus_name][1],
            "gen" => Dict{Any,Any}(),
            "load" => Dict{Any,Any}())
        bus_name_to_id[bus_names[i]] = i
    end

    p_data["bus"] = buses
    p_data["bus_name_to_id"] = bus_name_to_id

    csv_file_path = balance_file
    ts_df = CSV.read(csv_file_path, DataFrame)

    # Define a function to replace missing with 0 and find the maximum
    max_or_zero = x -> maximum(coalesce.(map(safe_parse_to_float, x), 0), init=0)

    # Group by 'Balancing Authority' and find the maximum 'Net Generation (MW) from Coal', assuming 0 for missing
    coal_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Coal") => max_or_zero => Symbol("Max Coal Generation (MW)"))
    ng_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Natural Gas") => max_or_zero => Symbol("Max Natural Gas Generation (MW)"))
    nuclear_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Nuclear") => max_or_zero => Symbol("Max Nuclear Generation (MW)"))
    hydro_ps_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Hydropower and Pumped Storage") => max_or_zero => Symbol("Max Hydro and Pumped Storage (MW)"))
    solar_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Solar") => max_or_zero => Symbol("Max Solar Generation (MW)"))
    wind_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Wind") => max_or_zero => Symbol("Max Wind Generation (MW)"))

    add_generators_from_df(coal_generation_per_authority, "coal", p_data)
    add_generators_from_df(ng_generation_per_authority, "natural_gas", p_data)
    add_generators_from_df(nuclear_generation_per_authority, "nuclear", p_data)
    add_generators_from_df(hydro_ps_generation_per_authority, "hydro", p_data)
    add_generators_from_df(solar_generation_per_authority, "solar", p_data)
    add_generators_from_df(wind_generation_per_authority, "wind", p_data)

    p_data["branch"] = Dict{Any,Any}()
    for (key, value) in branches
        p_data["branch"][value["branch_id"]] = value

        f_bus_id = bus_name_to_id[p_data["branch"][value["branch_id"]]["f_bus"]]
        t_bus_id = bus_name_to_id[p_data["branch"][value["branch_id"]]["t_bus"]]

        lat1 = p_data["bus"][f_bus_id]["lat"]
        lon1 = p_data["bus"][f_bus_id]["lon"]
        lat2 = p_data["bus"][t_bus_id]["lat"]
        lon2 = p_data["bus"][t_bus_id]["lon"]

        
        p_data["branch"][value["branch_id"]]["branch_name"] = key
        p_data["branch"][value["branch_id"]]["f_bus"] = f_bus_id
        p_data["branch"][value["branch_id"]]["t_bus"] = t_bus_id
        p_data["branch"][value["branch_id"]]["distance"] = haversine_distance(lat1, lon1, lat2, lon2)
    end

    p_data["bus_name_to_id"] = bus_name_to_id

    # initialize arcs_from dictionary
    p_data["arcs_from"] = Dict{Any,Any}()

    for (key, value) in p_data["branch"]
        i = value["f_bus"]
        j = value["t_bus"]
        # add to the arcs_from dictionary
        if haskey(p_data["arcs_from"], i)
            push!(p_data["arcs_from"][i], key)
        else
            p_data["arcs_from"][i] = [key]
        end
        if haskey(p_data["arcs_from"], j)
            push!(p_data["arcs_from"][j], key)
        else
            p_data["arcs_from"][j] = [key]
        end
    end

    json_data = JSON.json(p_data)
    open("data/topology/n1/power_system_data.json", "w") do file
        write(file, json_data)
    end

    return p_data
end