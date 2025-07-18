using DataStructures

function benders_ptdf_write_to_csv(simdir, master_obj, theta_val, phi_val, y_val; lambda_val=nothing)
    filename = joinpath(simdir, "output", "benders_progress.csv")

    # Decompose y_val
    gamma_val, s_power_val, s_energy_val = y_val

    # Prepare the base data row
    data_dict = OrderedDict(
        "master_objective" => master_obj,
        "theta_val" => theta_val,
        "phi_val" => phi_val,
        "total_line_upgrades" => count(x -> x >=0.015, gamma_val),
        "total_storage_power" => sum(s_power_val),
        "total_storage_energy" => sum(s_energy_val),
        "total_storage_count" => count(x -> x >= 0.015, s_energy_val)
    )
    
    # Add lambda if provided
    if lambda_val !== nothing
        data_dict["lambda"] = lambda_val
    end
    
    # Create DataFrame from dictionary
    df = DataFrame([data_dict])
    
    # Check if file exists to determine if we need headers
    file_exists = isfile(filename)
    
    if file_exists
        # Read existing headers to ensure compatibility
        existing_headers = names(CSV.read(filename, DataFrame; limit=1))
        
         # Check if lambda column exists in the file but not in our data
        if "lambda" in existing_headers && !("lambda" in keys(data_dict))
            # Add empty lambda to maintain column structure
            df.lambda = [missing]
        end
        
        # Write to CSV, appending to it
        CSV.write(filename, df; append=true, header=false)
    else
        # New file, write with headers
        CSV.write(filename, df)
    end
end