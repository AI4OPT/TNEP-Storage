using CSV, DataFrames

function compare_energy_csvs(csv1, csv2)
    df1 = CSV.read(csv1, DataFrame)
    df2 = CSV.read(csv2, DataFrame)
    
    differences = []
    
    for col in names(df1)
        # Skip identifier columns
        if col ∈ ["Node_Index", "Node_Name", "Lat", "Lon", "Hour"]
            continue
        end
        
        # Find rows where values differ beyond tolerance
        for i in 1:nrow(df1)
            if !isapprox(df1[i,col], df2[i,col], atol=1e-6)
                push!(differences, (
                    node=df1[i,"Node_Index"],
                    hour=df1[i,"Hour"],
                    column=col,
                    value1=df1[i,col],
                    value2=df2[i,col],
                    diff=abs(df1[i,col] - df2[i,col])
                ))
            end
        end
    end
    
    if isempty(differences)
        println("No differences found between the CSV files")
        return DataFrame(node=Int[], hour=Int[], column=String[], 
                         value1=Float64[], value2=Float64[], diff=Float64[])
    else
        return DataFrame(differences) |> 
               df -> sort(df, :diff, rev=true)  # Sort by largest differences
    end
 end

function compare_branch_csvs(csv1, csv2)
    df1 = CSV.read(csv1, DataFrame)
    df2 = CSV.read(csv2, DataFrame)

    differences = []

    for col in names(df1)
        # Skip identifier columns
        if col ∈ ["Branch_Index", "Lat1", "Lon1", "Lat2", "Lon2", "Hour"]
            continue
        end

        # Find rows where values differ beyond tolerance
        for i in 1:nrow(df1)
            if !isapprox(df1[i, col], df2[i, col], atol=1e-6)
                push!(differences, (
                    branch=df1[i, "Branch_Index"],
                    hour=df1[i, "Hour"],
                    column=col,
                    value1=df1[i, col],
                    value2=df2[i, col],
                    diff=abs(df1[i, col] - df2[i, col])
                ))
            end
        end
    end

    if isempty(differences)
        return DataFrame(branch=Int[], hour=Int[], column=String[], value1=Float64[], value2=Float64[], diff=Float64[])
    end

    # Convert to DataFrame, sort by largest differences, and return top 20
    return DataFrame(differences) |>
           df -> sort(df, :diff, rev=true)
end
