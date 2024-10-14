using CSV
using DataFrames
using Plots

function get_storage_candidates(data, filepath)
    # read the energy csv and locate the lat/lon pairs that either have energy imbalance or renewable curtailment
    df = CSV.read(filepath, DataFrame)
    df = select(df, Not([:Node_Index, :Node_Name, :Hour]))
    grouped_df = combine(groupby(df, [:Lat, :Lon]), names(df, Not([:Lat, :Lon])) .=> sum)

    renewable_types = data["param"]["renewable_types"]

    # Calculate the renewable curtailment
    grouped_df[!, :Renewable_Curtailment] = sum(grouped_df[!, Symbol(renewable * "_production_sum")] for renewable in renewable_types) .-
    sum(grouped_df[!, Symbol(renewable * "_sum")] for renewable in renewable_types)

    # Select only the relevant columns to display
    final_df = select(grouped_df, [:Lat, :Lon, :Energy_Imbalance_sum, :Renewable_Curtailment])

    # Filter for lat/lon pairs that have energy imbalance or renewable curtailment
    final_df = filter(row -> row.Energy_Imbalance_sum < 0 || row.Renewable_Curtailment > 0, final_df)

    keys_in_final_df = Set{String}()
    for (key, value) in data["bus"]
        lat = value["lat"]
        lon = value["lon"]
        if any(row -> row.Lat == lat && row.Lon == lon, eachrow(final_df))
            push!(keys_in_final_df, key)
        end
    end

    return keys_in_final_df
end

function intersect_storage_candidates_original(data, no_upgrades_dir)
    final_cand = Set(keys(data["bus"]))
    for date in data["param"]["dates"]
        filepath = joinpath(no_upgrades_dir, "output", date, "energy.csv")
        cand = get_storage_candidates(data, filepath)
        final_cand = intersect(final_cand, cand)
    end

    return final_cand
end

function intersect_storage_candidates(data, no_upgrades_dir)
    if !haskey(data["param"], "candidate_count_threshold")
        return intersect_storage_candidates_original(data, no_upgrades_dir)
    end

    final_cand = Set()
    threshold = data["param"]["candidate_count_threshold"]
    bus_counts = storage_candidate_counts(data, no_upgrades_dir)

    for (bus, count) in bus_counts
        if count >= threshold
            push!(final_cand, bus)
        end
    end
    return final_cand
end

function storage_candidate_counts(data, no_upgrades_dir)
    bus_counts = Dict()
    for i in 1:length(data["param"]["dates"])
        date = data["param"]["dates"][i]
        prob = data["param"]["representative_prob"][i]
        filepath = joinpath(no_upgrades_dir, "output", date, "energy.csv")
        cand = get_storage_candidates(data, filepath)

        for bus in cand
            bus_counts[bus] = get(bus_counts, bus, 0) + prob
        end
    end
    return bus_counts
end

function print_bus_counts_thresholds(bus_counts)
    for threshold in 5:-1:1
        count = sum(value >= threshold for value in values(bus_counts))
        println("Counts >= $threshold: $count buses meet the threshold.")
    end
end

function plot_bus_counts(no_upgrades_dir)
    data = JSON.parsefile(joinpath(no_upgrades_dir, "data.json"))
    bus_counts = storage_candidate_counts(data, no_upgrades_dir)
    vals = collect(values(bus_counts))
    # Generate x-axis values (0 to 1, with a step of 0.01)
    x_vals = 0:0.01:0.99
    # Calculate y-axis values (number of keys where the value > x)
    y_vals = [count(v -> v > x, vals) for x in x_vals]
    plot(x_vals, y_vals, label="Number of busses with value > x", xlabel="Threshold value (0-1)", ylabel="Number of busses", legend=:topright)
    savefig(joinpath(no_upgrades_dir, "visual", "bus_counts_plot.png"))
end