using TOML
using JSON
using CSV
using DataFrames

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
            data["gen"]["$idx"]["profile"]["$rep_index"] = df[:, gen_id] * renewable_scale
        end
    end
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
            zone_id = string(data["bus"]["$i"]["zone_id"])
            # zone power demand (MW) is disaggregated to buses proportional to Pd
            demand = fdemand_df[:, zone_id] * (data["bus"]["$i"]["Pd"] / data["zone_pd"][zone_id])
            data["bus"]["$i"]["load"]["$rep_index"] = demand

            # -- Populate the wind of all busses --
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

    # Calculate and check the sum of representatives' probabilities
    sum_probabilities = sum(values(data["representative_prob"]))
    if abs(sum_probabilities - 1) > 1e-5
        error_message = "Representatives' probabilities do not sum to 1: $sum_probabilities"
        throw(InvalidProbabilitySumError(error_message))
    end
    return data
end
