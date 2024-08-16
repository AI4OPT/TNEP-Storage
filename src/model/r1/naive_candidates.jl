using CSV
using DataFrames

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
    keys_not_in_final_df = Set{String}()
    for (key, value) in data["bus"]
        lat = value["lat"]
        lon = value["lon"]
        if any(row -> row.Lat == lat && row.Lon == lon, eachrow(final_df))
            push!(keys_in_final_df, key)
        else
            push!(keys_not_in_final_df, key)
        end
    end

    return keys_in_final_df, keys_not_in_final_df
end
