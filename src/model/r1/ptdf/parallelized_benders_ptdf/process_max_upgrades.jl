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

# This function will get the initial core point for the benders (an over-invested first stage)
function compute_superset_core_point(superdir)
    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    storage_energy_size = get(toml_data, "storage_energy_size", 2.5)

    # Get the initial core points (solved without benders)
    # has rep day optimal investments
    initial_optima_dir = toml_data["initial_optima_dir"]

    dates = toml_data["dates"]

    csv_trans_files = [joinpath(initial_optima_dir, a_date, "output", "line_investments.csv") for a_date in dates]
    csv_stor_files = [joinpath(initial_optima_dir, a_date, "output", "storage_investments.csv") for a_date in dates]

    trans_df = process_csv_max_upgrade(csv_trans_files, joinpath(superdir, "line_investments.csv"))
    stor_df = process_csv_max_storage(csv_stor_files, joinpath(superdir, "storage_investments.csv"))

    gamma_val = trans_df[:, :Upgrade_Lvl]
    s_energy_int_val = stor_df[:, :Storage_Energy] * 1 / storage_energy_size
    return ceil.(Int, gamma_val), ceil.(Int, s_energy_int_val)
end