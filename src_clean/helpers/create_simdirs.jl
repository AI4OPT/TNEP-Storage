using TOML

function create_simdirs(superdir::String, toml_data::Dict; only_feasibility=nothing)
    simdirs = String[]
    years = get(toml_data, "years", [toml_data["decarbonization_year"]])
    
    for rep in toml_data["dates"]
        for year in years
            simdir = joinpath(superdir, string(year) * rep[5:end])
            
            # Create directory if needed and check if it was created
            dir_existed = isdir(simdir)
            mkpath(simdir)
            
            # Create modified config
            config = deepcopy(toml_data)
            config["dates"] = [rep]
            config["representative_prob"] = [1.0]
            config["num_representatives"] = 1
            config["decarbonization_year"] = year
            
            if !isnothing(only_feasibility)
                config["only_feasibility"] = only_feasibility
            end
            
            # Write config
            open(joinpath(simdir, "config.toml"), "w") do io
                TOML.print(io, config)
            end
            
            push!(simdirs, simdir)
            
            # Only print if directory was actually created
            if !dir_existed
                println("Created simulation directory: $simdir")
            end
        end
    end
    
    return simdirs
end

function clear_directory(dirpath::String)
    """Clear all files in a directory if it exists."""
    if isdir(dirpath)
        for file in readdir(dirpath, join=true)
            if isfile(file)
                rm(file)
            end
        end
    end
end