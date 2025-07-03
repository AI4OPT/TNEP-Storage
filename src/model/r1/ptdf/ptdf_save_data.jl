using JuMP
using CSV, DataFrames

function save_tracked_constraints(simdir, model::JuMP.Model, n_violated)
    # Convert dictionary to DataFrame
    df = DataFrame(
        arc = [k[1] for k in keys(model.ext[:tracked_constraints])],
        rep = [k[2] for k in keys(model.ext[:tracked_constraints])],
        time = [k[3] for k in keys(model.ext[:tracked_constraints])],
        ub = [k[4] for k in keys(model.ext[:tracked_constraints])],
        tracked = collect(values(model.ext[:tracked_constraints])),
        violations_left = [n_violated for i in keys(model.ext[:tracked_constraints])]
    )

    # Sort DataFrame by arc first, then rep, then time
    # sort!(df, [:arc, :rep, :time])

    # Save to CSV
    CSV.write(joinpath(simdir, "output", "tracked_constraints.csv"), df)
end

function save_solve_time(simdir, solve_time)
    # Create the DataFrame with time information
    df = DataFrame(
        Variable = ["time_to_solve"],
        Value = [solve_time]
    )
    # Ensure output directory exists
    output_dir = joinpath(simdir, "output")
    mkpath(output_dir)
    # Save to CSV
    CSV.write(joinpath(output_dir, "time.csv"), df)
end

function save_power_injections(simdir, model, data)
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
            push!(rows, (rep=r, time=t, bus=bus, gen=g, variable="pg", value=pg_values[r, g, t]))
        end
        
        # Unserved energy
        for i in 1:N
            if ue_values[r, i, t] > 1e-6
                push!(rows, (rep=r, time=t, bus=i, gen=0, variable="ue", value=ue_values[r, i, t]))
            end
        end
        
        # Storage
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
    
    # Save to CSV for each representative period
    for r in 1:R
        datestring = data["param"]["dates"][r]
        output_dir = joinpath(simdir, "output", datestring)
        mkpath(output_dir)
        
        df_rep = filter(row -> row.rep == r, df)
        CSV.write(joinpath(simdir, "output", "power_injections.csv"), df_rep)
    end
end