using TOML

function create_simdirs(superdir::String, toml_data::Dict; only_feasibility=nothing)
    simdirs = String[]
    year = toml_data["decarbonization_year"]
    
    for rep in toml_data["dates"]
        simdir = joinpath(superdir, string(year) * rep[5:end])
        
        # Create directory if needed
        mkpath(simdir)
        
        # Create modified config
        config = deepcopy(toml_data)
        config["dates"] = [rep]
        config["representative_prob"] = [1.0]
        config["num_representatives"] = 1

        if !isnothing(only_feasibility)
            config["only_feasibility"] = only_feasibility
        end
        
        # Write config
        open(joinpath(simdir, "config.toml"), "w") do io
            TOML.print(io, config)
        end
        
        push!(simdirs, simdir)
        println("Created simulation directory: $simdir")
    end
    
    return simdirs
end