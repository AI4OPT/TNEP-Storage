using CSV

function process_csv_max_upgrade(csv_files::Vector{String}, output_file::String)
    """
    Process multiple CSV files and create a new CSV with maximum upgrade levels.
    """
    
    # Read all CSV files and concatenate
    all_dfs = [CSV.read(file, DataFrame) for file in csv_files]
    combined_df = vcat(all_dfs...)
    
    # Group by Branch_Index and take row with maximum upgrade level
    result_df = combine(groupby(combined_df, :Branch_Index)) do group
        max_idx = argmax(group.Upgrade_Lvl)
        return group[max_idx, :]
    end
    
    # Sort and save
    sort!(result_df, :Branch_Index)
    # CSV.write(output_file, result_df)
    
    return result_df
end

function process_csv_max_storage(csv_files::Vector{String}, output_file::String)
    """
    Process multiple CSV files and create a new CSV with maximum storage energy.
    """
    
    # Read all CSV files and concatenate
    all_dfs = [CSV.read(file, DataFrame) for file in csv_files]
    combined_df = vcat(all_dfs...)
    
    # Group by Node_Index and take maximum energy for each node
    result_df = combine(groupby(combined_df, :Node_Index)) do group
        max_idx = argmax(group.Storage_Energy)          
        return group[max_idx, :]
    end
    
    # Sort and save
    sort!(result_df, :Node_Index)
    # CSV.write(output_file, result_df)
    
    return result_df
end

function compute_superset_core_point(superdir; is_multistage::Bool=false)
    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    
    if is_multistage
        # Multistage: return Dict{Int => Tuple{Vector, Vector}}
        years = get(toml_data, "years", [toml_data["decarbonization_year"]])
        initial_optima_dirs = toml_data["initial_optima_dir"]  # Array of directories
        
        @assert length(initial_optima_dirs) == length(years) "Number of initial_optima_dirs must match number of years"
        
        result = Dict{Int, Vector{Vector{Float64}}}()
        
        for (year, initial_optima_dir) in zip(years, initial_optima_dirs)
            og_dates = toml_data["dates"]
            dates = [string(year) * a_date[5:end] for a_date in og_dates]
            csv_trans_files = [joinpath(initial_optima_dir, a_date, "output", "line_investments.csv") for a_date in dates]
            csv_stor_files = [joinpath(initial_optima_dir, a_date, "output", "storage_investments.csv") for a_date in dates]
            
            trans_df = process_csv_max_upgrade(csv_trans_files, joinpath(superdir, "line_investments_$(year).csv"))
            stor_df = process_csv_max_storage(csv_stor_files, joinpath(superdir, "storage_investments_$(year).csv"))
            
            gamma_val = trans_df[:, :Upgrade_Lvl]
            s_energy_val = stor_df[:, :Storage_Energy]
            
            gamma_val = safe_ceil.(gamma_val)
            s_energy_val = safe_ceil.(s_energy_val)
            
            # Check for manual override for this year
            if haskey(toml_data, "inv_dir")
                inv_dir = toml_data["inv_dir"]
                # Try year-specific file first, then fall back to generic
                trans_file = joinpath(inv_dir, "line_investments_$(year).csv")
                stor_file = joinpath(inv_dir, "storage_investments_$(year).csv")
                
                if isfile(trans_file) && isfile(stor_file)
                    trans_df = CSV.read(trans_file, DataFrame)
                    stor_df = CSV.read(stor_file, DataFrame)
                    gamma_val = trans_df[:, :Upgrade_Lvl]
                    s_energy_val = stor_df[:, :Storage_Energy]
                end
            end
            
            result[year] = [gamma_val, s_energy_val]
        end
        
        # Enforce monotonicity: later years >= earlier years
        sorted_years = sort(collect(keys(result)))
        for i in 2:length(sorted_years)
            prev_year = sorted_years[i-1]
            curr_year = sorted_years[i]
            
            # Enforce monotonicity for transmission investments (gamma)
            prev_gamma = result[prev_year][1]
            curr_gamma = result[curr_year][1]
            result[curr_year][1] = max.(curr_gamma, prev_gamma)
            
            # Enforce monotonicity for storage investments (s_energy)
            prev_storage = result[prev_year][2]
            curr_storage = result[curr_year][2]
            result[curr_year][2] = max.(curr_storage, prev_storage)
        end
        
        return result
    else
        # Single-stage: return Tuple{Vector, Vector} (existing logic)
        initial_optima_dir = toml_data["initial_optima_dir"]  # Single directory
        og_dates = toml_data["dates"]
        year = toml_data["decarbonization_year"]
        dates = [string(year) * a_date[5:end] for a_date in og_dates]
        
        csv_trans_files = [joinpath(initial_optima_dir, a_date, "output", "line_investments.csv") for a_date in dates]
        csv_stor_files = [joinpath(initial_optima_dir, a_date, "output", "storage_investments.csv") for a_date in dates]
        
        trans_df = process_csv_max_upgrade(csv_trans_files, joinpath(superdir, "line_investments.csv"))
        stor_df = process_csv_max_storage(csv_stor_files, joinpath(superdir, "storage_investments.csv"))
        
        gamma_val = trans_df[:, :Upgrade_Lvl]
        s_energy_val = stor_df[:, :Storage_Energy]
        gamma_val = safe_ceil.(gamma_val)
        s_energy_val = safe_ceil.(s_energy_val)
        
        if haskey(toml_data, "inv_dir")
            inv_dir = toml_data["inv_dir"]
            trans_file = joinpath(inv_dir, "line_investments.csv")
            stor_file = joinpath(inv_dir, "storage_investments.csv")
            trans_df = CSV.read(trans_file, DataFrame)
            stor_df = CSV.read(stor_file, DataFrame)
            gamma_val = trans_df[:, :Upgrade_Lvl]
            s_energy_val = stor_df[:, :Storage_Energy]
        end
        
        return (gamma_val, s_energy_val)
    end
end

function safe_ceil(x, tolerance=1e-8)
    # If x is very close to an integer, round to that integer
    # Otherwise, take the ceiling
    rounded = round(Int, x)
    if abs(x - rounded) < tolerance
        return Int(rounded)
    else
        return Int(ceil(x))
    end
end

function compare_y_vals(y_val_1, y_val_2; atol=1e-10, rtol=1e-8)
    gamma1 = y_val_1[1]
    gamma2 = y_val_2[1]
    s_energy_1 = y_val_1[2]
    s_energy_2 = y_val_2[2]
    
    # Use approximate comparison for floating point values
    gamma_diff = sum(abs.(gamma1 .- gamma2))
    energy_diff = sum(abs.(s_energy_1 .- s_energy_2))
    
    total_diff = gamma_diff + energy_diff
    
    # Return 0 if difference is within tolerance, otherwise return the difference
    return total_diff < atol + rtol * max(
        sum(abs.(gamma1)) + sum(abs.(s_energy_1)),
        sum(abs.(gamma2)) + sum(abs.(s_energy_2))
    ) ? 0.0 : total_diff
end