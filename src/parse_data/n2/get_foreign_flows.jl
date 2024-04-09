using JSON
using CSV
using DataFrames
using Missings
using Plots

INTERCHANGE_FILE = "data/topology/n1/EIA930_INTERCHANGE_2016_Jan_Jun.csv"
const EIA_TAMU_FILE = "data/topology/n1/eia_tamu.json"

function get_foreign_flows(interchange_file, date)
    df = CSV.read(interchange_file, DataFrame)
    eia_tamu = JSON.parsefile(EIA_TAMU_FILE)

    foreign_set = Set()
    for i in keys(eia_tamu)
        if eia_tamu[i] == "Canada" || eia_tamu[i] == "Mexico"
            push!(foreign_set, i)
        end
    end

    filtered_df = filter(row -> row["Directly Interconnected Balancing Authority"] ∈ foreign_set, df)
    filtered_df = filter(row -> row["Data Date"] == date, filtered_df)

    foreign_flows = Dict{Tuple{String, String}, Vector{Int64}}()

    for row in eachrow(filtered_df)
        # Extract key components
        ba = row["Balancing Authority"]
        dia = row["Directly Interconnected Balancing Authority"]
        key = (dia, ba)

        # Handle missing interchange numbers by pushing 0
        interchange = ismissing(row["Interchange (MW)"]) ? 0 : parse(Int64, replace(row["Interchange (MW)"], "," => ""))
        interchange *= -1.0

        # Check if the key exists in the dictionary
        if haskey(foreign_flows, key)
            push!(foreign_flows[key], interchange)  # If key exists, append the interchange number
        else
            foreign_flows[key] = [interchange]  # If key does not exist, create a new array with the interchange number
        end
    end

    return foreign_flows
end

function plot_foreign_flows(interchange_file, date)
    foreign_flows = get_foreign_flows(interchange_file, date)
    # Determine global y-axis limits
    all_values = vcat(values(foreign_flows)...)
    global_min = minimum(all_values)
    global_max = maximum(all_values)

    # Generate and save a plot for each key
    for (key, values) in foreign_flows
        p = plot(values, label="Interchange", legend=:topright, title="From $(key[1]) to $(key[2])", 
                ylims=(global_min, global_max), xlabel="Time", ylabel="MW")
        savefig(p, "plot_$(key[1])_$(key[2]).png")
    end
end