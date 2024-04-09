using JSON
using CSV
using DataFrames
include("../../helpers/haversine_distance.jl")

const DATA_FILE = "data/topology/tamu/power_system_data.json"
const GEOS_FILE = "data/geojson/world.geojson"

# Function to calculate the center of a MultiPolygon
function calculate_center(coordinates)
    total_x, total_y, count = 0.0, 0.0, 0
    for polygon in coordinates
        for ring in polygon
            for coord in ring
                total_x += coord[1]
                total_y += coord[2]
                count += 1
            end
        end
    end
    return (total_x / count, total_y / count)
end

# Function to read GeoJSON and return a dictionary of zoneName to center
function geojson_to_zone_centers(geojson_str)
    geojson = JSON.parsefile(geojson_str)
    centers = Dict()
    
    for feature in geojson["features"]
        zone_name = feature["properties"]["zoneName"]
        coordinates = feature["geometry"]["coordinates"]
        center = calculate_center(coordinates)
        centers[zone_name] = center
    end
    
    return centers
end

function merge_concat_dicts(dict1::Dict{String, Any}, dict2::Dict{String, Any})
    result = copy(dict1)  # Start with a copy of the first dictionary

    for (key, value) in dict2
        if haskey(result, key)
            # Concatenate the arrays if the key exists in both dictionaries
            result[key] = vcat(result[key], value)
        else
            # Simply copy the key-value pair if the key is unique to dict2
            result[key] = value
        end
    end

    return result
end

function aggregate_tamu()
    data = JSON.parsefile(DATA_FILE)

    centers = geojson_to_zone_centers(GEOS_FILE)

    aggregate = Dict(
        "per_unit" => false,
        "baseMVA" => 100.0, 
        "multinetwork" => false,
        "name" => "ISO_transport_model",
        "description" => "TAMU aggregated by ISO",
        "source_type" => "NREL Breakthrough Zenodo",
        "source_version" => "2021_Feb Version 0.4.2"
    )

    # aggregate busses first
    busses = Dict()
    iso2idx = Dict()
    bus_idx = 1
    for idx in keys(data["bus"])
        bus = data["bus"][idx]
        # create iso/bus entries
        if !haskey(iso2idx, bus["iso"])
            iso2idx[bus["iso"]] = bus_idx
            busses[bus_idx] = Dict(
                "bus_name" => bus["iso"],
                "lat" => centers[bus["iso"]][2],
                "lon" => centers[bus["iso"]][1],
                "zone_pd" => Dict(),
                "load" => Dict(),
                "gen" => Dict{String, Any}()
            )
            bus_idx += 1
        end
        # create zone and pd mappings for each iso/bus
        if !haskey(busses[iso2idx[bus["iso"]]]["zone_pd"], bus["zone_id"])
            busses[iso2idx[bus["iso"]]]["zone_pd"][bus["zone_id"]] = bus["Pd"]
        else
            busses[iso2idx[bus["iso"]]]["zone_pd"][bus["zone_id"]] += bus["Pd"]
        end
        # add generators
        busses[iso2idx[bus["iso"]]]["gen"] = merge_concat_dicts(busses[iso2idx[bus["iso"]]]["gen"], bus["gen"])
    end

    # create branches dictionary
    branches = Dict()
    branch2idx = Dict()
    branch_idx = 1
    inter_iso_df = DataFrame(idx = String[], lat1 = Float64[], lon1 = Float64[], lat2 = Float64[], lon2 = Float64[], inter_iso = Vector{Any}())
    for idx in keys(data["branch"])
        the_branch = data["branch"][idx]
        f_bus = the_branch["f_bus"]
        t_bus = the_branch["t_bus"]
        f_iso = data["bus"]["$f_bus"]["iso"]
        t_iso = data["bus"]["$t_bus"]["iso"]

        # if inter iso
        if f_iso != t_iso
            key = tuple(sort([f_iso, t_iso]))

            if !haskey(branch2idx, key) # create new branch entry
                branch2idx[key] = branch_idx
                f_iso_idx = iso2idx[f_iso]
                t_iso_idx = iso2idx[t_iso]
                dist = haversine_distance(busses[f_iso_idx]["lat"], busses[f_iso_idx]["lon"], busses[t_iso_idx]["lat"], busses[t_iso_idx]["lon"])

                branches[branch_idx] = Dict(
                    "f_bus" => f_iso_idx,
                    "t_bus" => t_iso_idx,
                    "rate_a" => the_branch["rate_a"],
                    "distance" => dist,
                    "num_branches" => 1
                )
                branch_idx += 1
            else # aggregate branch flow rates
                branches[branch2idx[key]]["rate_a"] = max(branches[branch2idx[key]]["rate_a"], the_branch["rate_a"])
                branches[branch2idx[key]]["num_branches"] += 1
            end
            # push to the inter_iso_df
            push!(inter_iso_df, (idx=idx, lat1=data["bus"]["$f_bus"]["lat"], lon1=data["bus"]["$f_bus"]["lon"], lat2=data["bus"]["$t_bus"]["lat"], lon2=data["bus"]["$t_bus"]["lon"], inter_iso=key))
        end
    end

    # makes arcs_from dictionary
    arcs_from = Dict()
    for idx in keys(branches)
        f_bus = branches[idx]["f_bus"]
        t_bus = branches[idx]["t_bus"]

        if !haskey(arcs_from, f_bus)
            arcs_from[f_bus] = [idx]
        else
            push!(arcs_from[f_bus], idx)
        end
        if !haskey(arcs_from, t_bus)
            arcs_from[t_bus] = [idx]
        else
            push!(arcs_from[t_bus], idx)
        end
    end

    # update gen dictionary
    for idx in keys(data["gen"])
        old_idx = data["gen"][idx]["gen_bus"]
        new_idx = iso2idx[data["bus"]["$old_idx"]["iso"]]
        data["gen"][idx]["gen_bus"] = new_idx
    end

    aggregate["bus"] = busses
    aggregate["branch"] = branches
    aggregate["arcs_from"] = arcs_from
    aggregate["zone_pd"] = data["zone_pd"]
    aggregate["idx2gen"] = data["idx2gen"]
    aggregate["gen"] = data["gen"]

    json_data = JSON.json(aggregate)
    open("data/topology/tamu/n2/tamu_aggregate_ps_maxes_data.json", "w") do file
        write(file, json_data)
    end

    CSV.write("data/topology/tamu/n2/inter_iso_branches.csv", inter_iso_df)

    return aggregate
end





