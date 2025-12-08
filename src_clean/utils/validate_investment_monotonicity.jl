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

function validate_prev(simdir)
    line_file_2 = joinpath(simdir, "output", "line_investments.csv")
    storage_file_2 = joinpath(simdir, "output", "storage_investments.csv")

    line_file_1 = joinpath(simdir, "previous_investment_dir", "line_investments.csv")
    storage_file_1 = joinpath(simdir, "previous_investment_dir", "storage_investments.csv")

    return compare_line_investment_growth(line_file_1, line_file_2) && compare_storage_investment_growth(storage_file_1, storage_file_2)

end