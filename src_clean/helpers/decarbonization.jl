using TOML
using JSON
using CSV
using DataFrames
using Dates

function update_decarbonization(simdir, data)
    # Load config file
    config_file = joinpath(simdir, "config.toml")
    toml_data = TOML.parsefile(config_file)

    d_data_file = toml_data["decarbonization"]
    year = Symbol(toml_data["decarbonization_year"])

    renewable_types = Set(toml_data["renewable_types"])
    nonrenewable_types = Set(toml_data["nonrenewable_types"])

    # Read decarbonization file
    decarbonization = CSV.read(d_data_file, DataFrame, delim=',')

    # Create the dictionary
    ratios = Dict()
    for row in eachrow(decarbonization)
        type_key = row[:Type]
        new_value = row[year]
        old_value = row[2]  # Assuming second column is the total count
        ratios[type_key] = new_value / old_value
    end

    # Scale up or down based off of generation type
    for (gen_id, gen) in data["gen"]
        gen_type = gen["gen_type"]

        if gen_type in keys(ratios)
            # scale up pmax and pmin
            gen["pmax"] *= ratios[gen_type]
            gen["pmin"] *= ratios[gen_type]
            if gen_type in renewable_types
                # scale up the time series profile too
                for (key, values) in gen["profile"]
                    gen["profile"][key] = map(x -> x * ratios[gen_type], values)  # Scale each element in the array
                end
            end
        end
    end

    # Now scale the load as well
    for (bus_id, bus) in data["bus"]
        for rep_index in keys(bus["load"])
            bus["load"][rep_index] *= ratios["load"]
        end
    end

    # # Scale up load
    # # TODO also scale up zone_pd, or delete it... it's not used anymore anyway
    # for (bus_id, bus) in data["bus"]
    #     for (key, values) in bus["load"]
    #         bus["load"][key] = map(x -> x * ratios["load"], values)

    return data
end


function edit_decarbonization()
    df = CSV.read("data/topology/tamu/decarbonization.csv", DataFrame)    

    # Find the row where the first column is "load"
    row_index = findfirst(df[!, 1] .== "load")

    # Modify the row starting from the second element
    if !isnothing(row_index)
        for j in 2:size(df, 2)
            df[row_index, j] = j == 2 ? df[row_index, j] : df[row_index, j - 1] * 1.03
        end
    end

    # Print the modified DataFrame to check
    println(df)

    # Optionally, save the modified DataFrame back to CSV
    CSV.write("data/topology/tamu/decarbonization.csv", df)
end
