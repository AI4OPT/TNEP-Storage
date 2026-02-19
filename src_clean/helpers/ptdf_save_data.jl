using JuMP
using CSV, DataFrames

function save_solve_time(simdir, solve_time)
    """
    Save solve time to CSV file, appending if file already exists.
    """
    # Ensure output directory exists
    filepath = joinpath(simdir, "output", "time.csv")

    # Create the new row
    df_new = DataFrame(
        Variable = ["time_to_solve"],
        Value = [solve_time]
    )

    # Append if file exists, otherwise write fresh
    if isfile(filepath)
        df_existing = CSV.read(filepath, DataFrame)
        df_combined = vcat(df_existing, df_new)
        CSV.write(filepath, df_combined)
    else
        CSV.write(filepath, df_new)
    end
end

function save_power_injections(simdir, model::JuMP.Model, data)
    """
    Save power injections (pg, ue, ch, dis) to CSV file.
    Works with both PTDFModel.jump_model and regular JuMP models.
    """
    # Get solution values
    pg_values = value.(model[:pg])
    ue_values = value.(model[:ue])
    ch_values = value.(model[:ch])
    
    # Check if dis variable exists and get its values
    dis_values = nothing
    has_dis = haskey(model, :dis)
    if has_dis
        dis_values = value.(model[:dis])
    end
    
    # Create storage for rows
    rows = []
    
    # Extract dimensions
    R = data["param"]["num_representatives"]
    G = length(data["gen"])
    N = length(data["bus"])
    T = data["param"]["num_hours"]
    
    # Gather all power injections
    for r in 1:R, t in 1:T
        # Generation
        for g in 1:G
            bus = data["gen"]["$g"]["gen_bus"]
            if abs(pg_values[r, g, t]) > 1e-6
                push!(rows, (rep=r, time=t, bus=bus, gen=g, variable="pg", value=pg_values[r, g, t]))
            end
        end
        
        # Unserved energy
        for i in 1:N
            if ue_values[r, i, t] > 1e-6
                push!(rows, (rep=r, time=t, bus=i, gen=0, variable="ue", value=ue_values[r, i, t]))
            end
        end
        
        # Charging
        for i in 1:N
            if abs(ch_values[r, i, t]) > 1e-6
                push!(rows, (rep=r, time=t, bus=i, gen=0, variable="ch", value=ch_values[r, i, t]))
            end
        end
        
        # Discharge (if dis variable exists)
        if has_dis
            for i in 1:N
                if abs(dis_values[r, i, t]) > 1e-6
                    push!(rows, (rep=r, time=t, bus=i, gen=0, variable="dis", value=dis_values[r, i, t]))
                end
            end
        end
    end
    
    # Convert to DataFrame
    df = DataFrame(rows)
    
    # Ensure output directory exists
    output_dir = joinpath(simdir, "output")
    mkpath(output_dir)
    
    # Save consolidated CSV
    CSV.write(joinpath(output_dir, "power_injections.csv"), df)
    
    # Optionally save per representative period if dates exist
    if haskey(data["param"], "dates")
        for r in 1:R
            datestring = data["param"]["dates"][r]
            date_output_dir = joinpath(output_dir, datestring)
            mkpath(date_output_dir)
            df_rep = filter(row -> row.rep == r, df)
            CSV.write(joinpath(date_output_dir, "power_injections.csv"), df_rep)
        end
    end
end