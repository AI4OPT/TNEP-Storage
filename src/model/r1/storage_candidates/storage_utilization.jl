function get_storage_metrics(full_df, storage_node)
    df = filter(row -> row.Node_Index == storage_node, full_df)

    # calculate average utilization
    total_throughput = sum(df.Charge .+ df.Discharge)
    average_utilization = total_throughput / nrow(df)

    # calculate time-based utilization
    active_hours = sum((df.Charge .> 0) .| (df.Discharge .> 0))
    time_based_utilization = active_hours / nrow(df)

    # Calculate peak utilization (maximum of charge or discharge observed)
    peak_utilization = maximum(df.Charge .+ df.Discharge)

    return [average_utilization, time_based_utilization, peak_utilization]
end

function analyze_storage_utilization(simdir)
    datafile = joinpath(simdir, "data.json")
    data = JSON.parsefile(datafile)

    storage_investments = CSV.read(joinpath(simdir, "output", "storage_investments.csv"), DataFrame)
    storage_nodes = storage_investments[storage_investments.Storage_Power .> 0, :Node_Index]

    storage_metrics = Dict(storage_node => [0.0, 0.0, 0.0] for storage_node in storage_nodes)

    # Iterate over each date with its representative probability
    for (i, each_date) in enumerate(data["param"]["dates"])
        prob = data["param"]["representative_prob"]["$i"]
        df = CSV.read(joinpath(simdir, "output", each_date, "energy.csv"), DataFrame)

        for storage_node in storage_nodes
            storage_metrics[storage_node] += get_storage_metrics(df, storage_node) * prob
        end
    end

    return storage_metrics
end

# sim/r1/weighted/rep_days/no_upgrades/2035

function assess_utilization_heuristic(simdir, no_upgrades_dir; start_ind=nothing, end_ind=nothing)
    datafile = joinpath(simdir, "data.json")
    data = JSON.parsefile(datafile)

    storage_investments = joinpath(simdir, "storage_investments.csv")
    storage_utilizations = Dict{Int, Vector{Float64}}()

    final_cand = Set()
    if start_ind !== nothing && end_ind !== nothing
        final_cand = Set(start_ind:end_ind)
    else
        final_cand = intersect_storage_candidates_original(data, no_upgrades_dir)
    end

    for cand in final_cand
        cand = parse(Int, string(cand))

        # consider storage only at the candidate location
        df = CSV.read(storage_investments, DataFrame)
        df[:, :Storage_Power] .= 0.0
        df[:, :Storage_Energy] .= 0.0
        df[cand, :Storage_Power] = 15.0
        df[cand, :Storage_Energy] = 15.0
        CSV.write(storage_investments, df)

        # run the operational model
        model, data = run_model(simdir)
        obj = objective_value(model)

        storage_metrics = analyze_storage_utilization(simdir)

        # Append each metric and the objective value to storage_utilizations
        for (node, metrics) in storage_metrics
            if haskey(storage_utilizations, node)
                append!(storage_utilizations[node], metrics)
                push!(storage_utilizations[node], obj)
            else
                storage_utilizations[node] = vcat(metrics, [obj])
            end
        end

        # Convert storage_utilizations to a DataFrame and save it
        utilizations_df = DataFrame(
            Node_Index = collect(keys(storage_utilizations)),
            Average_Utilization = [storage_utilizations[node][1] for node in keys(storage_utilizations)],
            Time_Based_Utilization = [storage_utilizations[node][2] for node in keys(storage_utilizations)],
            Peak_Utilization = [storage_utilizations[node][3] for node in keys(storage_utilizations)],
            Objective = [storage_utilizations[node][4] for node in keys(storage_utilizations)]
        )
        CSV.write(joinpath(simdir, "output", "storage_utilizations_candidates.csv"), utilizations_df)
    end

    return storage_utilizations
end

function process_csv(filepath::String)
    # Read the CSV file into a DataFrame
    df = CSV.read(filepath, DataFrame)
    
    # Sort the DataFrame by the Objective column in ascending order
    sort!(df, :Objective)
    
    # Round all Float64 columns to 3 decimal places
    for col in names(df)
        if eltype(df[!, col]) <: Float64
            df[!, col] = round.(df[!, col], digits=3)
        end
    end
    
    # Write the modified DataFrame back to the original CSV file
    CSV.write(filepath, df)
end





