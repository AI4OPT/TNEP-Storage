using TOML

function convert_to_floats(col_data, date)
    # Create an empty array for the floats
    float_data = Float64[]

    for item in col_data
        if ismissing(item)
            # Handle the missing value; you can also choose to append a default value
            println("Missing value encountered in " * date)
            push!(float_data, 0.0)
        elseif item isa AbstractString
            # Convert the string to a float and append to the float_data array
            item = replace(item, "," => "")
            push!(float_data, parse(Float64, item))
        else
            push!(float_data, Float64(item))
        end
    end
    return float_data
end

# this function should not alter the data["param"] field whatsoever!!!
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

    rep_index = 1

    for date in data["param"]["dates"]
        fdemand_df = demand_df[occursin.(date, string.(demand_df[!, "UTC Time at End of Hour"])), :]
        fsolar_df = solar_df[occursin.(date, string.(solar_df[!, "UTC Time at End of Hour"])), :]
        fwind_df = wind_df[occursin.(date, string.(wind_df[!, "UTC Time at End of Hour"])), :]

        # Initialize load for all busses to be 0
        for i in 1:length(data["bus"])
            data["bus"]["$i"]["load"]["$rep_index"] = zeros(num_h)
        end

        # Populate the load of all busses
        for i in 2:length(names(demand_df))
            col_name = names(demand_df)[i]
            col_data = fdemand_df[1:num_h, col_name]
            float_data = convert_to_floats(col_data, date)
            bus_id = data["bus_name_to_id"][col_name]
            data["bus"]["$bus_id"]["load"]["$rep_index"] = float_data * load_scale
        end

        # Populate the solar/wind profiles of all generators
        for i in 1:length(data["bus"])
            if "solar" in keys(data["bus"]["$i"]["gen"])
                gen_id = data["bus"]["$i"]["gen"]["solar"][1]
                bus_id = data["gen"]["$gen_id"]["gen_bus"]
                bus_name = data["bus"]["$bus_id"]["bus_name"]
                col_data = fsolar_df[1:num_h, bus_name]
                float_data = convert_to_floats(col_data, date)
                data["gen"]["$gen_id"]["profile"]["$rep_index"] = float_data * renewable_scale
            end

            if "wind" in keys(data["bus"]["$i"]["gen"])
                gen_id = data["bus"]["$i"]["gen"]["wind"][1]
                bus_id = data["gen"]["$gen_id"]["gen_bus"]
                bus_name = data["bus"]["$bus_id"]["bus_name"]
                col_data = fwind_df[1:num_h, bus_name]
                float_data = convert_to_floats(col_data, date)
                data["gen"]["$gen_id"]["profile"]["$rep_index"] = float_data * renewable_scale
            end
        end
        rep_index += 1
    end

     # Calculate and check the sum of representatives' probabilities
     sum_probabilities = sum(values(data["representative_prob"]))
     if abs(sum_probabilities - 1) > 1e-5
         error_message = "Representatives' probabilities do not sum to 1: $sum_probabilities"
         throw(InvalidProbabilitySumError(error_message))
     end

    return data
end

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

function add_params(simdir::String, data::Dict)
    # Load config file
    config_file = joinpath(simdir, "config.toml")
    toml_data = TOML.parsefile(config_file)

    # Add config parameters to data file
    data["param"] = toml_data

    return data
end

function add_params_profiles(simdir)
    data = add_params(simdir)

    # Add profiles to the data file
    data = add_profiles(data)

    # Save the updated data file to the simdir
    json_data = JSON.json(data)
    open(joinpath(simdir, "ps_profiled_data.json"), "w") do file
        write(file, json_data)
    end
    return data
end