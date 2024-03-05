# this function will load in the demand, wind, and solar time-series

using Dates

const v1_ts_dir = "data/time_series/v1/"

function create_time_series(balance_file, output_dir)
    csv_file_path = balance_file
    ts_df = CSV.read(csv_file_path, DataFrame)
    # Replace missing values in 'Demand (MW)' with values from 'Demand (MW) (Imputed)'
    ts_df[!, Symbol("Demand (MW)")] = coalesce.(ts_df[!, Symbol("Demand (MW)")], ts_df[!, Symbol("Demand (MW) (Imputed)")], 0.0)
    ts_df[!, Symbol("Net Generation (MW) from Solar")] = coalesce.(ts_df[!, Symbol("Net Generation (MW) from Solar")], ts_df[!, Symbol("Net Generation (MW) from Solar (Imputed)")], 0.0)
    ts_df[!, Symbol("Net Generation (MW) from Wind")] = coalesce.(ts_df[!, Symbol("Net Generation (MW) from Wind")], ts_df[!, Symbol("Net Generation (MW) from Wind (Imputed)")], 0.0)

    date_format = "mm/dd/yyyy I:MM:SS p" # Define the date format
    ts_df[!, Symbol("UTC Time at End of Hour")] = [DateTime(x, date_format) for x in ts_df[!, Symbol("UTC Time at End of Hour")]]

    demand_df = unstack(ts_df, Symbol("UTC Time at End of Hour"), Symbol("Balancing Authority"), Symbol("Demand (MW)"))
    solar_df = unstack(ts_df, Symbol("UTC Time at End of Hour"), Symbol("Balancing Authority"), Symbol("Net Generation (MW) from Solar"))
    wind_df = unstack(ts_df, Symbol("UTC Time at End of Hour"), Symbol("Balancing Authority"), Symbol("Net Generation (MW) from Wind"))

    sort!(demand_df, Symbol("UTC Time at End of Hour"))
    sort!(solar_df, Symbol("UTC Time at End of Hour"))
    sort!(wind_df, Symbol("UTC Time at End of Hour"))

    demand_csv_path = joinpath(output_dir, "demand.csv") 
    solar_csv_path = joinpath(output_dir, "solar.csv") 
    wind_csv_path = joinpath(output_dir, "wind.csv") 

    # Write the DataFrame to a CSV file
    CSV.write(demand_csv_path, demand_df)
    CSV.write(solar_csv_path, solar_df)
    CSV.write(wind_csv_path, wind_df)
end






