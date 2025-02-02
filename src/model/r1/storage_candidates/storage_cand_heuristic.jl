using LIBSVM
using DataFrames, CSV, Statistics

function get_node_shed(simdir)
    datafile = joinpath(simdir, "data.json")
    data = JSON.parsefile(datafile)

    # Initialize an empty dictionary for cumulative results
    cumulative_summary_dict = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64}}()

    # Iterate over each date with its representative probability
    for (i, each_date) in enumerate(data["param"]["dates"])
        prob = data["param"]["representative_prob"]["$i"]
        df = CSV.read(joinpath(simdir, "output", each_date, "energy.csv"), DataFrame)
        
        # Calculate wind and solar curtailment
        df[:, :wind_curtailment] = df[:, :wind_production] .- df[:, :wind]
        df[:, :solar_curtailment] = df[:, :solar_production] .- df[:, :solar]

        # Group by (Lat, Lon) and sum each column
        df = select(df, Not([:Node_Index, :Node_Name, :Hour]))
        grouped_df = combine(groupby(df, [:Lat, :Lon]), names(df, Not([:Lat, :Lon])) .=> sum)

        # Create node_summary_dict for this date, weighted by `prob`
        node_summary_dict = Dict(
            (row.Lat, row.Lon) => (
                row.Energy_Imbalance_sum * prob,
                row.wind_curtailment_sum * prob,
                row.solar_curtailment_sum * prob
            ) for row in eachrow(grouped_df)
        )
        
        # Update cumulative_summary_dict with weighted sums
        for (key, value) in node_summary_dict
            if haskey(cumulative_summary_dict, key)
                # If the key already exists, sum the values
                cumulative_summary_dict[key] = (
                    cumulative_summary_dict[key][1] + value[1],
                    cumulative_summary_dict[key][2] + value[2],
                    cumulative_summary_dict[key][3] + value[3]
                )
            else
                # If the key doesn't exist, initialize it
                cumulative_summary_dict[key] = value
            end
        end
    end

    rounded_dict = Dict(key => (round.(value, digits=3)) for (key, value) in cumulative_summary_dict)
    return rounded_dict
end

function get_storage_inv(simdir)
    df = CSV.read(joinpath(simdir, "output", "storage_investments.csv"), DataFrame)
    node_location_dict = Dict(row.Node_Index => (row.Lat, row.Lon) for row in eachrow(df))
    nodes_with_storage = df[(df.Storage_Power .> 0) .| (df.Storage_Energy .> 0), :]

    lat_lon_inv_set = Set()

    for node in nodes_with_storage[:, :Node_Index]
        push!(lat_lon_inv_set, node_location_dict[node])
    end

    return lat_lon_inv_set
end

function analyze_storage_inv(no_upgrades_dir, transmission_dir, inv_dir)
    no_upgrades_info = get_node_shed(no_upgrades_dir)
    transmission_info = get_node_shed(transmission_dir)
    lat_lon_inv_set = get_storage_inv(inv_dir)

    for i in lat_lon_inv_set
        println(no_upgrades_info[i], transmission_info[i])
    end

    for i in keys(no_upgrades_info)
        if i ∉ lat_lon_inv_set
            if no_upgrades_info[i] !== (0.0, 0.0, 0.0)
                println(no_upgrades_info[i], transmission_info[i])
            end
        end
    end

    # Prepare DataFrame with shed_info values and an investment indicator
    data = DataFrame(
        Lat = [key[1] for key in keys(shed_info)],
        Lon = [key[2] for key in keys(shed_info)],
        shed_val1 = [value[1] for value in values(shed_info)],
        shed_val2 = [value[2] for value in values(shed_info)],
        shed_val3 = [value[3] for value in values(shed_info)],
        investment = [key in lat_lon_inv_set for key in keys(shed_info)]  # True if in investment set, else False
    )

    # Convert `investment` to an integer (1 for True, 0 for False)
    data.investment = Int.(data.investment)

    # Prepare data for SVM
    X = hcat(data.shed_val1, data.shed_val2, data.shed_val3)'
    y = data.investment

    # Fit the SVM model with a linear kernel
    svm_model = svmtrain(X, y)
    
    # Test model on the other half of the data.
    (y_pred, decision_values) = svmpredict(svm_model, X)

    # Compute accuracy
    println(mean((y_pred .== y))*100)


    # Evaluate SVM separability (print model coefficients for interpretability)
    println("SVM Model Trained.")
    println("SVM Support Vectors:")
    println(svm_model.sv)
end

