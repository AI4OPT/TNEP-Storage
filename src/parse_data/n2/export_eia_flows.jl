using JSON
using CSV
using DataFrames
using Missings

INTERCHANGE_FILE = "EIA/EIA930_INTERCHANGE_2016_Jan_Jun.csv"
INTERCHANGE_FILES = [
    "EIA/EIA930_INTERCHANGE_2016_Jan_Jun.csv",
    "EIA/EIA930_INTERCHANGE_2016_Jul_Dec.csv"
]
THE_DATE = "01/12/2016"
PS_DATA = "data/topology/tamu/n2/tamu_aggregate_ps_data.json"

function resolve_eia_tamu(eia_tamu, name)
    ans = eia_tamu[name]

    if ans == "Canada" || ans == "Mexico"
        return nothing
    elseif startswith(ans, "inside_")
        return ans[8:end]
    elseif startswith(ans, "retired_")
        return ans[9:end]
    else
        return ans
    end
end

# create mapping of ISO name to index
function get_name2idx(data)
    name2idx = Dict()
    for i in keys(data["bus"])
        name2idx[data["bus"][i]["bus_name"]] = i
    end
    return name2idx
end


function export_eia_flows(INTERCHANGE_FILE, THE_DATE, PS_DATA)
    df = CSV.read(INTERCHANGE_FILE, DataFrame)
    eia_tamu = JSON.parsefile("data/topology/n1/eia_tamu.json")
    eia_coords = JSON.parsefile("data/topology/n1/xycoords.json")
    data = JSON.parsefile(PS_DATA)

    name2idx = get_name2idx(data)

    # filter for the specified date
    filtered_df = filter(row -> row["Data Date"] == THE_DATE, df)

    # create flows dictionary
    flows = Dict()

    for row in eachrow(filtered_df)
        from = row["Balancing Authority"]
        to = row["Directly Interconnected Balancing Authority"]
        hour = row["Hour Number"]
        amount = ismissing(row["Interchange (MW)"]) ? "0" : replace(row["Interchange (MW)"], "," => "")
        amount = parse(Int, amount)
        
        from = resolve_eia_tamu(eia_tamu, from)
        to = resolve_eia_tamu(eia_tamu, to)

        if isnothing(from) || isnothing(to)
            continue
        end

        # Check alphabetical order and adjust amount if necessary
        if from > to
            amount *= -1
            from, to = to, from  # Swap to ensure from and to are in alphabetical order
        end

        if !haskey(flows, (from, to))
            flows[(from,to)] = zeros(24)
        end
        flows[(from,to)][hour] += amount
    end

    flows_df = DataFrame(
        Lat1 = Float64[],
        Lon1 = Float64[],
        Lat2 = Float64[],
        Lon2 = Float64[],
        Hour = Int[],
        Power_Flow = Float64[]
    )

    for i in keys(flows)
        f_bus = name2idx[i[1]]
        t_bus = name2idx[i[2]]
        lat1 = data["bus"][f_bus]["lat"]
        lon1 = data["bus"][f_bus]["lon"]
        lat2 = data["bus"][t_bus]["lat"]
        lon2 = data["bus"][t_bus]["lon"]

        for j in 1:24
            push!(flows_df, (
                Lat1 = lat1,
                Lon1 = lon1,
                Lat2 = lat2,
                Lon2 = lon2,
                Hour = j,
                Power_Flow = flows[i][j] * 0.5
            ))
        end
    end

    the_date = replace(THE_DATE, "/" => "_")
    CSV.write("data/eia_cross_validation/flows_$the_date.csv", flows_df)

end



function max_abs_flow_per_pair(interchange_file, ps_data)
    df = CSV.read(interchange_file, DataFrame)
    select!(df, Not(["Local Time at End of Hour", "Region", "DIBA_Region"]))
    eia_tamu = JSON.parsefile("data/topology/n1/eia_tamu.json")
    data = JSON.parsefile(ps_data)

    name2idx = get_name2idx(data)

    df[!, "Directly Interconnected Balancing Authority"] = String.(df[:, "Directly Interconnected Balancing Authority"])
    df[!, "Balancing Authority"] = String.(df[:, "Balancing Authority"])
    df[!, "Interchange (MW)"] = replace.(coalesce.(df[!, "Interchange (MW)"], "0"), "," => "")
    df[!, "Interchange (MW)"] = parse.(Int, df[!, "Interchange (MW)"])


    for row in eachrow(df)
        from = row["Balancing Authority"]
        to = row["Directly Interconnected Balancing Authority"]
        hour = row["Hour Number"]

        from = resolve_eia_tamu(eia_tamu, from)
        to = resolve_eia_tamu(eia_tamu, to)

        if isnothing(from) || isnothing(to)
            continue
        end

        # replace with the tamu names
        row["Balancing Authority"] = from
        row["Directly Interconnected Balancing Authority"] = to
    end

    gdf = groupby(df, ["Balancing Authority", "Directly Interconnected Balancing Authority", "UTC Time at End of Hour"])
    sum_df = combine(gdf, "Interchange (MW)" => sum => "Total Interchange (MW)")

    return df
end

function set_max_flows(df, flows)
    for row in eachrow(df)
        from = row["Balancing Authority"]
        to = row["Directly Interconnected Balancing Authority"]
        amount = row["Interchange (MW)"]

        if from > to
            from, to = to, from  # Swap to ensure from and to are in alphabetical order
        end

        if !haskey(flows, (from, to))
            flows[(from,to)] = 0
        end
        flows[(from,to)] = max(amount, flows[(from,to)])
    end

    return flows
end

function get_max_flows(interchange_files, ps_data)
    flows = Dict()
    for filename in interchange_files
        df = max_abs_flow_per_pair(filename, ps_data)
        flows = set_max_flows(df, flows)
    end
    return flows
end

function set_eia_maxes(ps_data, flows)
    data = JSON.parsefile(ps_data)

    changes = Dict()

    for i in keys(data["branch"])
        f_bus = data["branch"][i]["f_bus"]
        t_bus = data["branch"][i]["t_bus"]
        from = data["bus"]["$f_bus"]["bus_name"]
        to = data["bus"]["$t_bus"]["bus_name"]
        pair = sort((from, to))
        if haskey(flows, pair)
            changes[pair] = (data["branch"][i]["rate_a"], flows[pair])
            data["branch"][i]["rate_a"] = flows[pair]
        end
    end

    json_data = JSON.json(data)
    open("data/topology/tamu/n2/tamu_aggregate_ps_eia_maxes_data.json", "w") do file
        write(file, json_data)
    end

    return changes
end