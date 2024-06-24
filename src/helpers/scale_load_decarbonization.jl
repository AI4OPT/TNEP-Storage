using CSV
using DataFrames

# Function to read, scale, and write CSV data
function scale_load_decarbonization(file_path::String, scale_factor::Float64)
    # Read the CSV file into a DataFrame
    df = CSV.read(file_path, DataFrame)

    # Identify the row where Type is 'load'
    load_row = findfirst(df.Type .== "load")

    # Scale the values in the row corresponding to 'load' for each year by the scale_factor
    initial_value = df[load_row, 2]  # The initial value to start scaling
    for col in 3:ncol(df)
        df[load_row, col] = initial_value * scale_factor^(col-2)
    end

    # Write the updated DataFrame back to the same CSV file
    CSV.write(file_path, df)
end