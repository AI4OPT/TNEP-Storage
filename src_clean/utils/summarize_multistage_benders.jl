using CSV, DataFrames, SHA

function summarize_multistage_trans(superdir)
    dir_path = joinpath(superdir, "benders_output")
    files = filter(f -> startswith(f, "line_investments_") && endswith(f, ".csv"), 
                   readdir(dir_path))
    
    # Parse iteration and year from filenames
    data = Dict()  # For sums
    config_data = Dict()  # For full configurations
    
    for file in files
        parts = split(replace(file, ".csv" => ""), "_")
        iter, year = parse(Int, parts[3]), parse(Int, parts[4])
        
        df = CSV.read(joinpath(dir_path, file), DataFrame)
        
        # Store sum
        data[(iter, year)] = sum(df.Upgrade_Lvl)
        
        # Store configuration as a sorted tuple of (Branch_Index, Upgrade_Lvl) for non-zero upgrades
        # This creates a unique fingerprint of the configuration
        upgrades = [(row.Branch_Index, row.Upgrade_Lvl) for row in eachrow(df) if row.Upgrade_Lvl > 0]
        sort!(upgrades)
        config_data[(iter, year)] = upgrades
    end
    
    # Build summary DataFrame
    iters = sort(unique([k[1] for k in keys(data)]))
    years = sort(unique([k[2] for k in keys(data)]))
    
    summary = DataFrame(Iteration = iters)
    
    # Add sum columns for each year
    for year in years
        summary[!, Symbol("$year")] = [get(data, (i, year), 0.0) for i in iters]
    end
    
    # Add configuration match columns for each year
    # These show if the configuration matches the previous iteration
    for year in years
        match_col = Symbol("$(year)_Match")
        matches = String[]
        
        for i in 1:length(iters)
            iter = iters[i]
            if i == 1
                push!(matches, "-")  # First iteration has nothing to compare to
            else
                prev_iter = iters[i-1]
                current_config = get(config_data, (iter, year), [])
                prev_config = get(config_data, (prev_iter, year), [])
                
                if current_config == prev_config
                    push!(matches, "✓")  # Exact match
                else
                    # Check if it's a subset or superset
                    if issubset(Set(current_config), Set(prev_config))
                        push!(matches, "⊂")  # Subset
                    elseif issubset(Set(prev_config), Set(current_config))
                        push!(matches, "⊃")  # Superset
                    else
                        push!(matches, "✗")  # Different
                    end
                end
            end
        end
        summary[!, match_col] = matches
    end
    
    return summary
end

function summarize_multistage_storage(superdir)
    dir_path = joinpath(superdir, "benders_output")
    files = filter(f -> startswith(f, "storage_investments_") && endswith(f, ".csv"), 
                   readdir(dir_path))
    
    # Parse iteration and year from filenames
    data = Dict()  # For sums
    config_data = Dict()  # For full configurations
    
    for file in files
        parts = split(replace(file, ".csv" => ""), "_")
        iter, year = parse(Int, parts[3]), parse(Int, parts[4])
        
        df = CSV.read(joinpath(dir_path, file), DataFrame)
        
        # Store sum of Storage_Energy
        data[(iter, year)] = sum(df.Storage_Energy)
        
        # Store configuration as a sorted tuple of (Node_Index, Storage_Energy) for non-zero storage
        storage_investments = [(row.Node_Index, row.Storage_Energy) 
                               for row in eachrow(df) if row.Storage_Energy > 0]
        sort!(storage_investments)
        config_data[(iter, year)] = storage_investments
    end
    
    # Build summary DataFrame
    iters = sort(unique([k[1] for k in keys(data)]))
    years = sort(unique([k[2] for k in keys(data)]))
    
    summary = DataFrame(Iteration = iters)
    
    # Add sum columns for each year
    for year in years
        summary[!, Symbol("$year")] = [get(data, (i, year), 0.0) for i in iters]
    end
    
    # Add configuration match columns for each year
    for year in years
        match_col = Symbol("$(year)_Match")
        matches = String[]
        
        for i in 1:length(iters)
            iter = iters[i]
            if i == 1
                push!(matches, "-")  # First iteration has nothing to compare to
            else
                prev_iter = iters[i-1]
                current_config = get(config_data, (iter, year), [])
                prev_config = get(config_data, (prev_iter, year), [])
                
                if current_config == prev_config
                    push!(matches, "✓")  # Exact match
                else
                    # Check if it's a subset or superset
                    if issubset(Set(current_config), Set(prev_config))
                        push!(matches, "⊂")  # Subset
                    elseif issubset(Set(prev_config), Set(current_config))
                        push!(matches, "⊃")  # Superset
                    else
                        push!(matches, "✗")  # Different
                    end
                end
            end
        end
        summary[!, match_col] = matches
    end
    
    return summary
end

# Combined function that returns both summaries
function summarize_multistage_investments(superdir)
    trans_summary = summarize_multistage_trans(superdir)
    storage_summary = summarize_multistage_storage(superdir)
    
    return (transmission = trans_summary, storage = storage_summary)
end

# Usage examples:
# trans_summary = summarize_multistage_trans("path/to/superdir")
# storage_summary = summarize_multistage_storage("path/to/superdir")
# 
# Or get both at once:
# summaries = summarize_multistage_investments("path/to/superdir")
# println(summaries.transmission)
# println(summaries.storage)