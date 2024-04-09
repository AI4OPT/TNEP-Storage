using TOML
using JSON
using CSV
using DataFrames
using Dates

include("../../parse_data/n2/get_foreign_flows.jl")

function add_params(simdir::String)
    # Load config file
    config_file = joinpath(simdir, "config.toml")
    toml_data = TOML.parsefile(config_file)

    # Load power systems data file
    p_data_file = toml_data["power_system_data"]
    data = JSON.parsefile(p_data_file)

    # Add config parameters to data file
    data["param"] = toml_data

    return data
end

function update_generation_profile(data, rep_index, bus_id, resource_type, df, renewable_scale)
    # Check if the resource type exists in the current bus's generators
    if resource_type in keys(data["bus"][bus_id]["gen"])
        # Get generator IDs for the specific resource type
        gens = data["bus"][bus_id]["gen"][resource_type]
        # Update the profile for each generator ID
        for idx in gens
            gen_id = string(data["idx2gen"]["$idx"])
            if !haskey(data["gen"]["$idx"], "profile")
                data["gen"]["$idx"]["profile"] = Dict()
            end
            data["gen"]["$idx"]["profile"]["$rep_index"] = df[:, gen_id] * renewable_scale
        end
    end
end

function update_load(data, rep_index, bus_id, df, load_scale)
    demand = zeros(data["param"]["num_hours"])
    for zone_id in keys(data["bus"][bus_id]["zone_pd"])
        demand += df[:, zone_id] * data["bus"][bus_id]["zone_pd"][zone_id] / data["zone_pd"][zone_id]
    end
    data["bus"][bus_id]["load"]["$rep_index"] = demand * load_scale
end

function add_profiles(data)
    num_r = data["param"]["num_representatives"]
    num_h = data["param"]["num_hours"]
    load_scale = data["param"]["load_scale"]
    renewable_scale = data["param"]["load_scale"]

    # Check if number of representatives is equal to the dates
    if num_r != length(data["param"]["dates"])
        throw(ArgumentError("Number of representatives does not match the config file"))
    end

    # Create the representative probabilities
    data["representative_prob"] = Dict{Int, Float64}()
    for i in 1:num_r
        data["representative_prob"][i] = 1 / num_r
    end

    # Load the demand, solar, and wind time-series
    ts_dir = data["param"]["time_series_dir"]
    demand_df = CSV.read(joinpath(ts_dir, "demand.csv"), DataFrame)
    solar_df = CSV.read(joinpath(ts_dir, "solar.csv"), DataFrame)
    wind_df = CSV.read(joinpath(ts_dir, "wind.csv"), DataFrame)
    hydro_df = CSV.read(joinpath(ts_dir, "hydro.csv"), DataFrame)

    rep_index = 1

    for date in data["param"]["dates"]
        fdemand_df = demand_df[occursin.(date, string.(demand_df[!, "UTC Time"])), :]
        fsolar_df = solar_df[occursin.(date, string.(solar_df[!, "UTC"])), :]
        fwind_df = wind_df[occursin.(date, string.(wind_df[!, "UTC"])), :]
        fhydro_df = hydro_df[occursin.(date, string.(hydro_df[!, "UTC"])), :]
        
        for i in keys(data["bus"])
            # -- Populate the load of all busses --
            update_load(data, rep_index, i, fdemand_df, load_scale)

            # -- Populate the wind of all busses --
            update_generation_profile(data, rep_index, i, "wind_offshore", fwind_df, renewable_scale)
            update_generation_profile(data, rep_index, i, "wind", fwind_df, renewable_scale)
            update_generation_profile(data, rep_index, i, "solar", fsolar_df, renewable_scale)
            update_generation_profile(data, rep_index, i, "hydro", fhydro_df, renewable_scale)
        end
        rep_index += 1
    end
end

function add_params_profiles(simdir)
    data = add_params(simdir)
    add_profiles(data)
    add_foreign_imports(data)

    # Calculate and check the sum of representatives' probabilities
    sum_probabilities = sum(values(data["representative_prob"]))
    if abs(sum_probabilities - 1) > 1e-5
        error_message = "Representatives' probabilities do not sum to 1: $sum_probabilities"
        throw(InvalidProbabilitySumError(error_message))
    end
    return data
end

function add_foreign_imports(data)
    if data["param"]["foreign_imports"] == false
        return
    end

    tamu2idx = Dict()
    for i in keys(data["bus"])
        tamu2idx[data["bus"][i]["bus_name"]] = i
    end

    eia_pair2idx = Dict()
    gen_idx = length(data["gen"]) + 1
    rep_idx = 1

    eia_tamu = JSON.parsefile(EIA_TAMU_FILE)
    for date in data["param"]["dates"]
        formatted_date = Dates.format(Date(date, "yyyy-mm-dd"), "mm/dd/yyyy")
        foreign_flows = get_foreign_flows(data["param"]["eia_interchange"], formatted_date)

        for i in keys(foreign_flows)
            if !haskey(eia_pair2idx, i)
                eia_pair2idx[i] = gen_idx

                # create new generator at that node
                # first update the bus dictionary
                bus_idx = tamu2idx[eia_tamu[i[2]]]
                haskey(data["bus"][bus_idx]["gen"], "foreign") ? push!(data["bus"][bus_idx]["gen"]["foreign"], gen_idx) : data["bus"][bus_idx]["gen"]["foreign"] = [gen_idx]
                # then update the gen dictionary
                data["gen"]["$gen_idx"] = Dict(
                    "profile" => Dict(),
                    "cost" => [0, 0, 0],
                    "gen_type" => "foreign",
                    "status" => 1,
                    "model" => 2,
                    "gen_bus" => parse(Int, bus_idx),
                    "ncost" => 3
                )
                gen_idx += 1
            end
            a_gen_idx = eia_pair2idx[i]
            data["gen"]["$a_gen_idx"]["profile"]["$rep_idx"] = foreign_flows[i]
        end
        rep_idx += 1

    end
end