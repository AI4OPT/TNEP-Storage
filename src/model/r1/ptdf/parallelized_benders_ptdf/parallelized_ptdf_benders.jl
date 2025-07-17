using JuMP
using Gurobi
using Dualization
using CSV, DataFrames
using DataStructures

function parallelized_ptdf_benders(superdir)

    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    dates = toml_data["dates"]

    # Compute the initial core point
    compute_superset_core_point(superdir)

    # Create master_output
    if !isdir(joinpath(superdir, "output"))
        mkdir(joinpath(superdir, "output"))
    end
    if !isdir(joinpath(superdir, "benders_output"))
        mkdir(joinpath(superdir, "benders_output"))
    end

    # Create the subproblem directories
    for a_date in dates
        if !isdir(joinpath(superdir, a_date))
            mkdir(joinpath(superdir, a_date))
            mkdir(joinpath(superdir, a_date, "output"))
        end
    end

    # Set up the data by reading the config
    data = set_up_data(superdir)

    # Create master problem
    master, y, theta = define_master_ptdf(data)

end

function compute_superset_core_point(superdir)
    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)

    # Get the initial core points (solved without benders)
    # has rep day optimal investments
    initial_optima_dir = toml_data["initial_optima_dir"]

    dates = toml_data["dates"]

    csv_trans_files = [joinpath(initial_optima_dir, a_date, "output", "line_investments.csv") for a_date in dates]
    csv_stor_files = [joinpath(initial_optima_dir, a_date, "output", "storage_investments.csv") for a_date in dates]

    process_csv_max_upgrade(csv_trans_files, joinpath(superdir, "line_investments.csv"))
    process_csv_max_storage(csv_stor_files, joinpath(superdir, "storage_investments.csv"))
end



    
    

