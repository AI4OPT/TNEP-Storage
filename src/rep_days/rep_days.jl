using TOML
using JSON
using CSV
using DataFrames
using Dates
using Clustering
using Distances
using Statistics

toml_data = TOML.parsefile(config_file)
data = JSON.parsefile(toml_data["power_system_data"])
ts_dir = data["param"]["time_series_dir"]
demand_df = CSV.read(joinpath(ts_dir, "demand.csv"), DataFrame)
solar_df = CSV.read(joinpath(ts_dir, "solar.csv"), DataFrame)
wind_df = CSV.read(joinpath(ts_dir, "wind.csv"), DataFrame)

# make the wind df
wind_df[!, :total_wind] = zeros(size(wind_df, 1))
for (idx, gen) in data["gen"]
    if gen["gen_type"] == "wind"
        gen_id = data["idx2gen"][idx]
        wind_df[!, :total_wind] .= wind_df[!, :total_wind] .+ wind_df[!, Symbol(gen_id)]
    end
end

wind_df = wind_df[:, ["UTC", "total_wind"]]
wind_df[!, :Year] = year.(DateTime.(wind_df[!, :UTC], "yyyy-mm-dd HH:MM:SS"))
wind_df = select(wind_df, :UTC, :Year, :total_wind)
rename!(wind_df, :UTC => :Timestamp)

# make the solar df
solar_df[!, :total_solar] = zeros(size(solar_df, 1))
for (idx, gen) in data["gen"]
    if gen["gen_type"] == "solar"
        gen_id = data["idx2gen"][idx]
        solar_df[!, :total_solar] .= solar_df[!, :total_solar] .+ solar_df[!, Symbol(gen_id)]
    end
end

solar_df = solar_df[:, ["UTC", "total_solar"]]
solar_df[!, :Year] = year.(DateTime.(solar_df[!, :UTC], "yyyy-mm-dd HH:MM:SS"))
solar_df = select(solar_df, :UTC, :Year, :total_solar)
rename!(solar_df, :UTC => :Timestamp)

# make demand df
zone_columns = keys(data["zone_pd"])
demand_df[!, :total_demand] = sum(demand_df[:, Symbol(k)] for k in zone_columns)
demand_df = demand_df[:, ["UTC Time", "total_demand"]]
demand_df[!, :Year] = year.(DateTime.(demand_df[Symbol("UTC Time")], "yyyy-mm-dd HH:MM:SS"))
demand_df = select(demand_df, Symbol("UTC Time"), :Year, :total_demand)
rename!(demand_df, Symbol("UTC Time") => :Timestamp)

dir = "src/rep_days"

CSV.write(joinpath(dir, "demand_df.csv"), demand_df)
CSV.write(joinpath(dir, "wind_df.csv"), wind_df)
CSV.write(joinpath(dir, "solar_df.csv"), solar_df)

ts_input_data = load_timeseries_data("$dir", years=[2016], T=24)
clust_res = run_clust(ts_input_data;method="kmedoids",representation="medoid",n_clust=5)

ts_clust_data = clust_res.clust_data

demands = ts_clust_data.data["demand_df-total_demand"]
demand_df = CSV.read(joinpath(dir, "demand_df.csv"), DataFrame)

# Step 1: Extract daily demand patterns from demand_df
num_days = Int(size(demand_df, 1) / 24)  # Number of full days in the year
daily_patterns = [demand_df[24*(i-1)+1:24*i, :total_demand] for i in 1:num_days]

# Step 2: Compare each day in demand_df to each column in demands
matches = []

for j in 1:5  # For each day in demands (5 days)
    best_match_day = -1
    best_distance = Inf
    for i in 1:num_days  # For each day in demand_df
        daily_demand = daily_patterns[i]
        # Compute Euclidean distance between the clustering results and daily demand
        dist = sum((demands[:, j] - daily_demand).^2)
        if dist < best_distance
            best_distance = dist
            best_match_day = i
        end
    end
    push!(matches, best_match_day)
end

# Step 3: Output the matching days (timestamps)
matching_timestamps = demand_df[[24*(i-1)+1 for i in matches], :Timestamp]
println("Matching Days:")
println(matching_timestamps)

# -- DEMAND ONLY --
# ["2016-04-18 00:00:00", 
# "2016-07-12 00:00:00", 
# "2016-01-19 00:00:00", 
# "2016-10-14 00:00:00", 
# "2016-09-02 00:00:00"]

# ---- MY OWN REP DAYS IMPLEMENTATION
dir = "src/rep_days"
wind_df = CSV.read(joinpath(dir, "wind_df.csv"), DataFrame)
solar_df = CSV.read(joinpath(dir, "solar_df.csv"), DataFrame)
demand_df = CSV.read(joinpath(dir, "demand_df.csv"), DataFrame)

time_series = hcat(wind_df[!, :total_wind], solar_df[!, :total_solar], demand_df[!, :total_demand])'

function normalize_time_series(time_series)
    return (time_series .- mean(time_series, dims=2)) ./ std(time_series, dims=2)
end

function split_into_days(normalized_data)
    num_days = Int(size(normalized_data, 2) / 24)
    return reshape(normalized_data, 3, 24, num_days)
end

function kmedoids_clustering_with_weights(daily_data, k)
    # Flatten the data for each day and cluster
    num_days = size(daily_data, 3)
    flattened_data = hcat([vec(daily_data[:, :, i]) for i in 1:num_days]...)

    # Calculate pairwise Euclidean distance between each day's data
    dist_matrix = pairwise(Euclidean(), flattened_data)

    # Perform k-medoids clustering
    result = kmedoids(dist_matrix, k)

    # Calculate weights: the number of points in each cluster
    weights = [count(i -> i == cluster, result.assignments) for cluster in 1:k]

    # Assuming the medoids refer to days, map medoid index to the first hour of the corresponding day
    medoid_day_indices = result.medoids
    # Multiply by 24 to get the index of the first hour of that day
    medoid_hour_indices = (medoid_day_indices .- 1) .* 24 .+ 1  # Convert day index to the first hour of the day
    # Extract the corresponding timestamps for the first hour of each medoid day
    medoid_dates = wind_df[medoid_hour_indices, :Timestamp]

    normalized_weights = weights ./ num_days
    medoid_dates_only = Date.(DateTime.(medoid_dates, "yyyy-mm-dd HH:MM:SS"))
    date_weight_dict = Dict(string(medoid_dates_only[i]) => normalized_weights[i] for i in 1:k)

    return result.assignments, date_weight_dict
end

function kmedoids_wrapper(time_series, k)
    normalized_data = normalize_time_series(time_series)
    daily_data = split_into_days(normalized_data)
    assignments, date_dict =  kmedoids_clustering_with_weights(daily_data, k)
    format_for_toml_output(date_dict)
end

function format_for_toml_output(dict)
    # Extract dates and probabilities from the dictionary
    dates = collect(keys(dict))
    probabilities = round.(collect(values(dict)), digits=4)
    
    # Format the output
    println("dates = ", dates)
    println("representative_prob = ", probabilities)
end

`Dict{String, Float64} with 5 entries:
  "2016-07-27" => 0.120219
  "2016-04-24" => 0.172131
  "2016-09-08" => 0.199454
  "2016-10-25" => 0.221311
  "2016-10-01" => 0.286885`

`
dates = ["2016-10-01", "2016-07-06", "2016-06-11", "2016-03-19", "2016-03-29", "2016-01-14", "2016-09-13", "2016-09-09", "2016-08-15", "2016-08-08"]
representative_prob = [0.2377, 0.0219, 0.0847, 0.0874, 0.0847, 0.1557, 0.0874, 0.0628, 0.112, 0.0656]
`