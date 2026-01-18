using CSV, DataFrames, SHA

function compare_y_vals(y_val_1, y_val_2; atol=1e-10)
    gamma1 = y_val_1[1]
    gamma2 = y_val_2[1]
    s_energy_1 = y_val_1[2]
    s_energy_2 = y_val_2[2]
    
    total_diff = sum(abs.(gamma1 .- gamma2)) + sum(abs.(s_energy_1 .- s_energy_2))
    
    return total_diff < atol ? 0.0 : total_diff
end

function compare_configs(config1, config2; atol=1e-10)
    # Check if same length
    if length(config1) != length(config2)
        return false
    end
    
    # Extract indices and values separately
    if length(config1) == 0
        return true  # Both empty configs are equal
    end
    
    # Get indices and values
    indices1 = [c[1] for c in config1]
    indices2 = [c[2] for c in config2]
    values1 = [c[2] for c in config1]
    values2 = [c[2] for c in config2]
    
    # Indices must match exactly
    if indices1 != indices2
        return false
    end
    
    # Use compare_y_vals logic for values
    # Create y_val format: ([values], [])  - we only care about first component
    y_val_1 = (values1, [])
    y_val_2 = (values2, [])
    
    diff = compare_y_vals(y_val_1, y_val_2; atol=atol)
    
    return diff == 0.0
end

function summarize_multistage_trans(superdir; atol=1e-10)
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
        upgrades = [(row.Branch_Index, row.Upgrade_Lvl) for row in eachrow(df) if row.Upgrade_Lvl > atol]
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
    for year in years
        match_col = Symbol("$(year)_Match")
        matches = String[]
        
        for i in 1:length(iters)
            iter = iters[i]
            if i == 1
                push!(matches, "-")
            else
                prev_iter = iters[i-1]
                current_config = get(config_data, (iter, year), [])
                prev_config = get(config_data, (prev_iter, year), [])
                
                if compare_configs(current_config, prev_config; atol=atol)
                    push!(matches, "✓")
                else
                    push!(matches, "✗")
                end
            end
        end
        summary[!, match_col] = matches
    end
    
    return summary
end

function summarize_multistage_storage(superdir; atol=1e-10)
    dir_path = joinpath(superdir, "benders_output")
    files = filter(f -> startswith(f, "storage_investments_") && endswith(f, ".csv"), 
                   readdir(dir_path))
    
    # Parse iteration and year from filenames
    data = Dict()
    config_data = Dict()
    
    for file in files
        parts = split(replace(file, ".csv" => ""), "_")
        iter, year = parse(Int, parts[3]), parse(Int, parts[4])
        
        df = CSV.read(joinpath(dir_path, file), DataFrame)
        
        # Store sum of Storage_Energy
        data[(iter, year)] = sum(df.Storage_Energy)
        
        # Store configuration, filtering out near-zero values
        storage_investments = [(row.Node_Index, row.Storage_Energy) 
                               for row in eachrow(df) if row.Storage_Energy > atol]
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
                push!(matches, "-")
            else
                prev_iter = iters[i-1]
                current_config = get(config_data, (iter, year), [])
                prev_config = get(config_data, (prev_iter, year), [])
                
                if compare_configs(current_config, prev_config; atol=atol)
                    push!(matches, "✓")
                else
                    push!(matches, "✗")
                end
            end
        end
        summary[!, match_col] = matches
    end
    
    return summary
end

# Combined function that returns both summaries
function summarize_multistage_investments(superdir; atol=1e-10)
    trans_summary = summarize_multistage_trans(superdir; atol=atol)
    storage_summary = summarize_multistage_storage(superdir; atol=atol)
    
    return (transmission = trans_summary, storage = storage_summary)
end