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
    CSV.write(output_file, result_df)
    
    return result_df
end

function process_csv_max_storage(csv_files::Vector{String}, output_file::String)
    """
    Process multiple CSV files and create a new CSV with maximum storage power and energy.
    """
    
    # Read all CSV files and concatenate
    all_dfs = [CSV.read(file, DataFrame) for file in csv_files]
    combined_df = vcat(all_dfs...)
    
    # Group by Node_Index and take maximum power and energy for each node
    result_df = combine(groupby(combined_df, :Node_Index)) do group
        # Find the row with maximum total storage (power + energy)
        max_idx = argmax(group.Storage_Power .+ group.Storage_Energy)
        base_row = group[max_idx, :]
        
        # But use actual maximum values for power and energy separately
        base_row.Storage_Power = maximum(group.Storage_Power)
        base_row.Storage_Energy = maximum(group.Storage_Energy)
        
        return base_row
    end
    
    # Sort and save
    sort!(result_df, :Node_Index)
    CSV.write(output_file, result_df)
    
    return result_df
end