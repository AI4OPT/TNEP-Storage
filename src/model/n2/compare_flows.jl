using CSV
using DataFrames

function sort_flows(filename)
    df = CSV.read(filename, DataFrame)

    # Round float columns to 3 decimal places
    for col in names(df)
        if eltype(df[!, col]) <: AbstractFloat
            df[!, col] = round.(df[!, col], digits=3)
        end
    end
    
    sorted_df = sort(df, [:Lat1, :Lon1])
    CSV.write(filename, sorted_df)
end