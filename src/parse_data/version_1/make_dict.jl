using JSON
using CSV
using DataFrames

function add_generators_from_df(df, gen_type, p_data)
    for row in eachrow(df)
        bus = row[Symbol("Balancing Authority")]
        capacity = row[2]  # Assuming the max capacity is the third column

        if capacity > 0
            bus_id = String(bus)  # Convert bus name to string if it isn't already
            gen_id = gen_type * "_" * bus_id  # Create a unique generator ID

            # Add the generator
            p_data["gen"][gen_id] = Dict(
                "gen_bus" => bus_id, 
                "pmax" => capacity,
                "gen_type" => gen_type
            )
        end
    end
end


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

csv_file_path = "data/version_1/EIA930_INTERCHANGE_2019_Jan_Jun.csv"
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
    if !ismissing(row["Interchange (MW)"])
        parsed_value = tryparse(Float64, row["Interchange (MW)"])
        if parsed_value !== nothing
            interchange_value = parsed_value
        end
    end
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
            branches[pair_key] = Dict("f_bus" => pair[1], "t_bus" => pair[2], "branch_id" => length(branches) + 1)
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
buses = Dict{String,Any}()
for bus in bus_names
    buses[bus] = Dict("bus_type" => 1) # Add other necessary bus attributes here
end

csv_file_path = "data/version_1/EIA930_BALANCE_2019_Jan_Jun.csv"
ts_df = CSV.read(csv_file_path, DataFrame)

# Define a function to replace missing with 0 and find the maximum
safe_parse_to_float(x) = ismissing(x) ? 0.0 : tryparse(Float64, x) === nothing ? 0.0 : parse(Float64, x)
max_or_zero = x -> maximum(coalesce.(map(safe_parse_to_float, x), 0), init=0)

# Group by 'Balancing Authority' and find the maximum 'Net Generation (MW) from Coal', assuming 0 for missing
coal_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Coal") => max_or_zero => Symbol("Max Coal Generation (MW)"))
ng_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Natural Gas") => max_or_zero => Symbol("Max Natural Gas Generation (MW)"))
nuclear_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Nuclear") => max_or_zero => Symbol("Max Nuclear Generation (MW)"))
hydro_ps_generation_per_authority = combine(groupby(ts_df, Symbol("Balancing Authority")), Symbol("Net Generation (MW) from Hydropower and Pumped Storage") => max_or_zero => Symbol("Max Hydro and Pumped Storage (MW)"))

add_generators_from_df(coal_generation_per_authority, "coal", p_data)
add_generators_from_df(ng_generation_per_authority, "natural gas", p_data)
add_generators_from_df(nuclear_generation_per_authority, "nuclear", p_data)
add_generators_from_df(hydro_ps_generation_per_authority, "hydro", p_data)

p_data["bus"] = buses
p_data["branch"] = branches

json_data = JSON.json(p_data)
open("power_system_data.json", "w") do file
    write(file, json_data)
end