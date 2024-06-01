using JSON
using CSV
using DataFrames
include("../../helpers/haversine_distance.jl")
include("../../helpers/safe_parse_to_float.jl")

const DATA_DIR = "tamu/base_grid/"
const BUS_FILE = "bus.csv"
const BRANCH_FILE = "branch.csv"
const BUS2SUB_FILE = "bus2sub.csv"
# const SUB_FILE = "tamu/base_grid/sub.csv"
const SUB_FILE = "data/geojson/upgraded_sub.csv"
const GEN_FILE = "plant.csv"
const DEMAND_FILE = "demand.csv"
const GENCOST_FILE = "gencost.csv"

function load_df_and_filter_region(filepath, region)
    df = CSV.read(filepath, DataFrame)
    if !isnothing(region)
        df_filtered = df[df.interconnect .== region, :]
        return df_filtered
    else
        return df
    end
end

function get_lat_lon_sub(bus2sub, sub_info, bus_id)
    sub = bus2sub[bus_id]
    ans = sub_info[sub]
    return ans["lat"], ans["lon"], ans["name"], ans["iso"]
end

function create_p_data(region=nothing)
    p_data = Dict(
        "per_unit" => true,
        "baseMVA" => 100.0, 
        "multinetwork" => false,
        "name" => "TAMU_DC_model",
        "description" => "TAMU",
        "source_type" => "NREL Breakthrough Zenodo",
        "source_version" => "2021_Feb Version 0.4.2"
    )

    df_bus = load_df_and_filter_region(joinpath(DATA_DIR, BUS_FILE), region)
    df_sub = load_df_and_filter_region(SUB_FILE, region)
    df_bus2sub = load_df_and_filter_region(joinpath(DATA_DIR, BUS2SUB_FILE), region)
    df_branch = load_df_and_filter_region(joinpath(DATA_DIR, BRANCH_FILE), region)
    df_gen = load_df_and_filter_region(joinpath(DATA_DIR, GEN_FILE), region)
    df_gencost = load_df_and_filter_region(joinpath(DATA_DIR, GENCOST_FILE), region)

    # create bus2sub and sub_info dictionary
    bus2sub = Dict(row[:bus_id] => row[:sub_id] for row in eachrow(df_bus2sub))
    sub_info = Dict(row[:sub_id] => Dict("name" => row[:name], "lat" => row[:lat], "lon" => row[:lon], "iso" => row[:zoneName]) for row in eachrow(df_sub))

    # create busses dictionary
    # create zonal demand disaggregation dict
    busses = Dict()
    zone_pd = Dict()
    bus2idx = Dict()
    bus_idx = 1
    for row in eachrow(df_bus)
        lat, lon, sub_name, iso_name = get_lat_lon_sub(bus2sub, sub_info, row["bus_id"])

        bus2idx[row["bus_id"]] = bus_idx
        busses[bus_idx] = Dict(
            "bus_type" => row["type"],
            "load" => Dict(),
            "gen" => Dict(),
            "Pd" => row["Pd"],
            "zone_id" => row["zone_id"],
            "interconnect" => row["interconnect"],
            "lat" => lat,
            "lon" => lon,
            "bus_name" => sub_name,
            "sub" => bus2sub[row["bus_id"]],
            "iso" => iso_name
        )
        zone_pd[row["zone_id"]] = get(zone_pd, row["zone_id"], 0) + row["Pd"]
        bus_idx += 1
    end

    # create branches dictionary
    branches = Dict()
    idx2branch = Dict()
    arcs_from = Dict{Any,Any}()
    branch_idx = 1
    for row in eachrow(df_branch)
        f_bus = bus2idx[row["from_bus_id"]]
        t_bus = bus2idx[row["to_bus_id"]]

        dist = haversine_distance(busses[f_bus]["lat"], busses[f_bus]["lon"], busses[t_bus]["lat"], busses[t_bus]["lon"])
        dist = maximum((0.25, dist)) # chose 0.25 km if the branch connects two busses at the same sub

        branch = Dict(
            "f_bus" => f_bus,
            "t_bus" => t_bus,
            "br_r" => row["r"],
            "br_x" => row["x"],
            "rate_a" => row["rateA"],
            "br_status" => row["status"],
            "br_type" => row["branch_device_type"],
            "angmin" => row["angmin"],
            "angmax" => row["angmax"],
            "distance" => dist
        )

        # update arcs_from dictionary
        if !haskey(arcs_from, f_bus)
            arcs_from[f_bus] = [branch_idx]
        else
            push!(arcs_from[f_bus], branch_idx)
        end
        if !haskey(arcs_from, t_bus)
            arcs_from[t_bus] = [branch_idx]
        else
            push!(arcs_from[t_bus], branch_idx)
        end

        # update branch dictionary
        branches[branch_idx] = branch

        branch_idx += 1
    end

    # create gencosts dictionary
    gencosts = Dict(row[:plant_id] => Dict("type" => row[:type], 
                "n" => row[:n], 
                "coeffs" => [row[:c0], row[:c1], row[:c2]]) for row in eachrow(df_gencost))

    # create gen dictionary
    gens = Dict()
    idx2gen = Dict()
    gen_idx = 1
    for row in eachrow(df_gen)

        plant = row["plant_id"]
        gen_bus = bus2idx[row["bus_id"]]
        idx2gen[gen_idx] = plant

        gens[gen_idx] = Dict(
            "profile" => Dict(),
            "cost" => gencosts[plant]["coeffs"],
            "model" => gencosts[plant]["type"],
            "ncost" => gencosts[plant]["n"],
            "gen_bus" => gen_bus,
            "gen_type" => row["type"],
            "pmax" => row["Pmax"],
            "pmin" => row["Pmin"],
            "status" => row["status"]
        )

        # update busses to track gens
        if !haskey(busses[gen_bus]["gen"], row["type"])
            busses[gen_bus]["gen"][row["type"]] = [gen_idx]
        else
            push!(busses[gen_bus]["gen"][row["type"]], gen_idx)
        end

        gen_idx += 1
    end

    p_data["bus"] = busses
    p_data["branch"] = branches
    p_data["arcs_from"] = arcs_from
    p_data["gen"] = gens
    # p_data["bus2sub"] = bus2sub
    p_data["sub_info"] = sub_info
    p_data["zone_pd"] = zone_pd
    p_data["idx2gen"] = idx2gen

    json_data = JSON.json(p_data)
    open("data/topology/tamu/power_system_data.json", "w") do file
        write(file, json_data)
    end

    return p_data
end

        

    
    


