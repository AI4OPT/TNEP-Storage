using TOML
using JSON
using CSV
using DataFrames
using Dates

FROM_UTC_TIMEZONE_SHIFT = Dict(
    "CDT" => -5
)

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

function update_generation_profile(data, rep_index, bus_id, resource_type, df)
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
            data["gen"]["$idx"]["profile"]["$rep_index"] = df[:, gen_id]
        end
    end
end

function update_load(data, rep_index, bus_id, df)
    zone_id = string(data["bus"][bus_id]["zone_id"])
    demand = df[:, zone_id] * data["bus"][bus_id]["Pd"] / data["zone_pd"][zone_id]
    data["bus"][bus_id]["load"]["$rep_index"] = demand
end

function add_profiles(data)
    num_r = data["param"]["num_representatives"]
    num_h = data["param"]["num_hours"]

    specific_hours_range = get(data["param"], "specific_hours_range", nothing)
    if specific_hours_range !== nothing
        data["param"]["num_hours"] = specific_hours_range[2] - specific_hours_range[1] + 1
        num_h = data["param"]["num_hours"]
    end

    # Check if number of representatives is equal to the dates
    if num_r != length(data["param"]["dates"])
        throw(ArgumentError("Number of representatives does not match the config file"))
    end

    # Create the representative probabilities
    if !haskey(data["param"], "representative_prob") || data["param"]["representative_prob"] == true                  
        data["param"]["representative_prob"] = Dict{Int, Float64}()
        for i in 1:num_r
            data["param"]["representative_prob"][i] = 1 / num_r
        end
    end

    # Load the demand, solar, and wind time-series
    ts_dir = data["param"]["time_series_dir"]
    demand_df = CSV.read(joinpath(ts_dir, "demand.csv"), DataFrame)
    solar_df = CSV.read(joinpath(ts_dir, "solar.csv"), DataFrame)
    wind_df = CSV.read(joinpath(ts_dir, "wind.csv"), DataFrame)
    hydro_df = CSV.read(joinpath(ts_dir, "hydro.csv"), DataFrame)

    # Shift timezone if specified
    timezone = "UTC"
    if haskey(data["param"], "timezone")
        timezone = data["param"]["timezone"] 
    end

    if timezone != "UTC" && haskey(FROM_UTC_TIMEZONE_SHIFT, timezone)
        shift_hours = FROM_UTC_TIMEZONE_SHIFT[timezone]

        # Function to apply the shift
        function apply_timezone_shift!(df, time_column::Symbol)
            df[!, time_column] = DateTime.(df[!, time_column], "yyyy-mm-dd HH:MM:SS") .+ Hour(shift_hours)
        end

        apply_timezone_shift!(demand_df, Symbol("UTC Time"))
        apply_timezone_shift!(solar_df, Symbol("UTC"))
        apply_timezone_shift!(wind_df, Symbol("UTC"))
        apply_timezone_shift!(hydro_df, Symbol("UTC"))
    end

    rep_index = 1

    for date in data["param"]["dates"]
        fdemand_df = demand_df[occursin.(date, string.(demand_df[!, "UTC Time"])), :]
        fsolar_df = solar_df[occursin.(date, string.(solar_df[!, "UTC"])), :]
        fwind_df = wind_df[occursin.(date, string.(wind_df[!, "UTC"])), :]
        fhydro_df = hydro_df[occursin.(date, string.(hydro_df[!, "UTC"])), :]

        if specific_hours_range !== nothing
            hr = specific_hours_range[1]:specific_hours_range[2]
            fdemand_df = fdemand_df[hr, :]
            fsolar_df  = fsolar_df[hr, :]
            fwind_df   = fwind_df[hr, :]
            fhydro_df  = fhydro_df[hr, :]
        end

        for i in keys(data["bus"])
            # -- Populate the load of all busses --
            update_load(data, rep_index, i, fdemand_df)

            # -- Populate the wind of all busses --
            update_generation_profile(data, rep_index, i, "wind_offshore", fwind_df)
            update_generation_profile(data, rep_index, i, "wind", fwind_df)
            update_generation_profile(data, rep_index, i, "solar", fsolar_df)
            update_generation_profile(data, rep_index, i, "hydro", fhydro_df)
        end
        rep_index += 1
    end
end

"""
function make_date_subdirs(simdir, data)
    for date in data["param"]["dates"]
        dirname = joinpath(simdir, "output", date)
        if !isdir(dirname)
            mkdir(dirname)
        end
    end
end
"""

function add_params_profiles(simdir)
    data = add_params(simdir)
    add_profiles(data)
    # make_date_subdirs(simdir, data)

    # Calculate and check the sum of representatives' probabilities
    sum_probabilities = sum(values(data["param"]["representative_prob"]))
    if abs(sum_probabilities - 1) > 1e-3
        error_message = "Representatives' probabilities do not sum to 1: $sum_probabilities"
        throw(InvalidProbabilitySumError(error_message))
    end
    return data
end

struct InvalidProbabilitySumError <: Exception
    message::String
end
