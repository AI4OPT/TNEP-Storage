using DataStructures

function benders_ptdf_write_to_csv(output_dir, y_val, master_obj, theta_val, phi_val, total_ue)
    filename = joinpath(output_dir, "output", "benders_progress.csv")

    # Decompose y_val
    gamma_val, s_energy_int_val = y_val

    # Prepare the base data row
    data_dict = OrderedDict(
        "master_obj_lb" => master_obj,
        "theta_val" => theta_val,
        "phi_val" => phi_val,
        "master_obj_ub" => master_obj + phi_val - theta_val,
        "total_ue" => total_ue,
        "total_line_upgrades" => count(x -> x >=0.015, gamma_val),
        "sum_line_upgrades" => sum(gamma_val),
        "total_storage_energy" => sum(s_energy_int_val),
        "total_storage_count" => count(x -> x >=0.015, s_energy_int_val)
    )

    # Create DataFrame from dictionary
    df = DataFrame([data_dict])
    
    # Check if file exists to determine if we need headers
    file_exists = isfile(filename)
    
    if file_exists
        # Read existing headers to ensure compatibility
        existing_headers = names(CSV.read(filename, DataFrame; limit=1))
        
        # Write to CSV, appending to it
        CSV.write(filename, df; append=true, header=false)
    else
        # New file, write with headers
        CSV.write(filename, df)
    end
end