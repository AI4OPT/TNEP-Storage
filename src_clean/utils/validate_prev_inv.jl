using CSV, DataFrames

function compare_line_investment_growth(
    line_file_1::String, 
    line_file_2::String;
    tolerance::Float64=1e-6
)
    """
    Compare line investment levels between two CSV files.
    Returns true if all upgrade levels in file_2 are >= file_1 (within tolerance).
    """
    
    # Read and sort by Branch_Index
    df1 = CSV.read(line_file_1, DataFrame)
    df2 = CSV.read(line_file_2, DataFrame)
    sort!(df1, :Branch_Index)
    sort!(df2, :Branch_Index)
    
    # Validate alignment
    if nrow(df1) != nrow(df2) || df1.Branch_Index != df2.Branch_Index
        return false
    end
    
    # Check all upgrade levels: file_2 >= file_1 (within tolerance)
    return all(df2.Upgrade_Lvl .>= df1.Upgrade_Lvl .- tolerance)
end

function compare_storage_investment_growth(
    storage_file_1::String,
    storage_file_2::String;
    tolerance::Float64=1e-6
)
    """
    Compare storage investment levels between two CSV files.
    Returns true if all storage levels in file_2 are >= file_1 (within tolerance).
    """
    
    # Read and sort by Node_Index
    df1 = CSV.read(storage_file_1, DataFrame)
    df2 = CSV.read(storage_file_2, DataFrame)
    sort!(df1, :Node_Index)
    sort!(df2, :Node_Index)
    
    # Validate alignment
    if nrow(df1) != nrow(df2) || df1.Node_Index != df2.Node_Index
        return false
    end
    
    # Check all storage levels: file_2 >= file_1 (within tolerance)
    return all(df2.Storage_Energy .>= df1.Storage_Energy .- tolerance)
end

function validate_prev_inv(simdir)
    line_file_2 = joinpath(simdir, "output", "line_investments.csv")
    storage_file_2 = joinpath(simdir, "output", "storage_investments.csv")

    line_file_1 = joinpath(simdir, "previous_investment_dir", "line_investments.csv")
    storage_file_1 = joinpath(simdir, "previous_investment_dir", "storage_investments.csv")

    return compare_line_investment_growth(line_file_1, line_file_2) && compare_storage_investment_growth(storage_file_1, storage_file_2)

end

function compare_line_investment_differences(
    line_file_1::String, 
    line_file_2::String;
    tolerance::Float64=1e-6,
    show_all::Bool=false  # If true, show all lines; if false, only show differences
)
    """
    Print differences in line investment levels between two CSV files.
    Useful for debugging trust region behavior and investment progression.
    """
    # Read and sort by Branch_Index
    df1 = CSV.read(line_file_1, DataFrame)
    df2 = CSV.read(line_file_2, DataFrame)
    sort!(df1, :Branch_Index)
    sort!(df2, :Branch_Index)
    
    # Validate alignment
    if nrow(df1) != nrow(df2)
        println("ERROR: Files have different number of rows!")
        println("  File 1: $(nrow(df1)) rows")
        println("  File 2: $(nrow(df2)) rows")
        return
    end
    
    if df1.Branch_Index != df2.Branch_Index
        println("ERROR: Files have different Branch_Index values!")
        return
    end
    
    # Print header
    println("\n" * "="^80)
    println("Investment Comparison:")
    println("  File 1: $line_file_1")
    println("  File 2: $line_file_2")
    println("="^80)
    
    # Calculate differences
    diff = df2.Upgrade_Lvl .- df1.Upgrade_Lvl
    
    # Count changes
    n_increased = count(diff .> tolerance)
    n_decreased = count(diff .< -tolerance)
    n_unchanged = count(abs.(diff) .<= tolerance)
    
    println("\nSummary:")
    println("  Increased:  $n_increased lines")
    println("  Decreased:  $n_decreased lines")
    println("  Unchanged:  $n_unchanged lines")
    println("  Total:      $(nrow(df1)) lines")
    
    # Print detailed differences
    if show_all
        println("\nAll Lines:")
        println("  Branch_Index | File1_Lvl | File2_Lvl | Difference | Status")
        println("  " * "-"^70)
        for i in 1:nrow(df1)
            status = if abs(diff[i]) <= tolerance
                "SAME"
            elseif diff[i] > 0
                "INCREASED"
            else
                "DECREASED"
            end
            println("  $(lpad(df1.Branch_Index[i], 12)) | $(lpad(round(df1.Upgrade_Lvl[i], digits=4), 9)) | $(lpad(round(df2.Upgrade_Lvl[i], digits=4), 9)) | $(lpad(round(diff[i], digits=4), 10)) | $status")
        end
    else
        # Only show differences
        diff_indices = findall(abs.(diff) .> tolerance)
        
        if isempty(diff_indices)
            println("\nNo differences found (within tolerance=$tolerance)")
        else
            println("\nDifferences ($(length(diff_indices)) lines):")
            println("  Branch_Index | File1_Lvl | File2_Lvl | Difference | Change")
            println("  " * "-"^70)
            for i in diff_indices
                change = diff[i] > 0 ? "↑ INCREASE" : "↓ DECREASE"
                println("  $(lpad(df1.Branch_Index[i], 12)) | $(lpad(round(df1.Upgrade_Lvl[i], digits=4), 9)) | $(lpad(round(df2.Upgrade_Lvl[i], digits=4), 9)) | $(lpad(round(diff[i], digits=4), 10)) | $change")
            end
        end
    end
    
    println("="^80 * "\n")
    
    # Return summary stats
    return (
        n_increased=n_increased,
        n_decreased=n_decreased,
        n_unchanged=n_unchanged,
        max_increase=isempty(diff) ? 0.0 : maximum(diff),
        max_decrease=isempty(diff) ? 0.0 : minimum(diff)
    )
end

function compare_storage_investment_differences(
    storage_file_1::String, 
    storage_file_2::String;
    tolerance::Float64=1e-6,
    show_all::Bool=false  # If true, show all nodes; if false, only show differences
)
    """
    Print differences in storage investment levels between two CSV files.
    Useful for debugging trust region behavior and storage investment progression.
    """
    # Read and sort by Node_Index
    df1 = CSV.read(storage_file_1, DataFrame)
    df2 = CSV.read(storage_file_2, DataFrame)
    sort!(df1, :Node_Index)
    sort!(df2, :Node_Index)
    
    # Validate alignment
    if nrow(df1) != nrow(df2)
        println("ERROR: Files have different number of rows!")
        println("  File 1: $(nrow(df1)) rows")
        println("  File 2: $(nrow(df2)) rows")
        return
    end
    
    if df1.Node_Index != df2.Node_Index
        println("ERROR: Files have different Node_Index values!")
        return
    end
    
    # Print header
    println("\n" * "="^80)
    println("Storage Investment Comparison:")
    println("  File 1: $storage_file_1")
    println("  File 2: $storage_file_2")
    println("="^80)
    
    # Calculate differences
    diff = df2.Storage_Energy .- df1.Storage_Energy
    
    # Count changes
    n_increased = count(diff .> tolerance)
    n_decreased = count(diff .< -tolerance)
    n_unchanged = count(abs.(diff) .<= tolerance)
    
    println("\nSummary:")
    println("  Increased:  $n_increased nodes")
    println("  Decreased:  $n_decreased nodes")
    println("  Unchanged:  $n_unchanged nodes")
    println("  Total:      $(nrow(df1)) nodes")
    
    # Print detailed differences
    if show_all
        println("\nAll Nodes:")
        println("  Node_Index | File1_Energy | File2_Energy | Difference | Status")
        println("  " * "-"^75)
        for i in 1:nrow(df1)
            status = if abs(diff[i]) <= tolerance
                "SAME"
            elseif diff[i] > 0
                "INCREASED"
            else
                "DECREASED"
            end
            println("  $(lpad(df1.Node_Index[i], 10)) | $(lpad(round(df1.Storage_Energy[i], digits=2), 12)) | $(lpad(round(df2.Storage_Energy[i], digits=2), 12)) | $(lpad(round(diff[i], digits=2), 10)) | $status")
        end
    else
        # Only show differences
        diff_indices = findall(abs.(diff) .> tolerance)
        
        if isempty(diff_indices)
            println("\nNo differences found (within tolerance=$tolerance)")
        else
            println("\nDifferences ($(length(diff_indices)) nodes):")
            println("  Node_Index | File1_Energy | File2_Energy | Difference | Change")
            println("  " * "-"^75)
            for i in diff_indices
                change = diff[i] > 0 ? "↑ INCREASE" : "↓ DECREASE"
                println("  $(lpad(df1.Node_Index[i], 10)) | $(lpad(round(df1.Storage_Energy[i], digits=2), 12)) | $(lpad(round(df2.Storage_Energy[i], digits=2), 12)) | $(lpad(round(diff[i], digits=2), 10)) | $change")
            end
        end
    end
    
    println("="^80 * "\n")
    
    # Return summary stats
    return (
        n_increased=n_increased,
        n_decreased=n_decreased,
        n_unchanged=n_unchanged,
        max_increase=isempty(diff) ? 0.0 : maximum(diff),
        max_decrease=isempty(diff) ? 0.0 : minimum(diff),
        total_energy_change=sum(diff)
    )
end