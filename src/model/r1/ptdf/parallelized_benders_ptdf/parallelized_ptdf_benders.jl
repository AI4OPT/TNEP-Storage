using JuMP
using Gurobi
using CSV, DataFrames
using Base.Threads
include("../../run_model.jl")
include("../../../../helpers/compute_gen_cost.jl")
include("../base_ptdf.jl")

include("process_max_upgrades.jl")
include("master_problem.jl")
include("subproblem.jl")
include("benders_logging.jl")
include("stabilization.jl")
include("process_max_upgrades.jl")


# This function sets up the superdirectory and the subproblem directories, and generates each data.json file
# It also copies the master config into subproblem specific configs
function set_up_directories(superdir)
    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    dates = toml_data["dates"]
    rep_prob = toml_data["representative_prob"]

    initial_optima_dir = toml_data["initial_optima_dir"]

    # Create master_output
    mkpath(joinpath(superdir, "output"))
    mkpath(joinpath(superdir, "benders_output"))

    for a_date in dates
        # Create subproblem directory
        simdir = joinpath(superdir, a_date)
        mkpath(simdir)
        mkpath(joinpath(simdir, "output"))

        # Make subproblem specific config
        config = deepcopy(toml_data)
        config["dates"] = [a_date]
        config["representative_prob"] = [1.0]
        config["num_representatives"] = 1

        open(joinpath(simdir, "config.toml"), "w") do io
            TOML.print(io, config)
        end

        # Warm-start subproblem with tracked constraints
        # tracked_file = joinpath(initial_optima_dir, a_date, "output", "tracked_constraints.csv")
        # cp(tracked_file, joinpath(simdir, "tracked_constraints.csv"), force=true)

        # Create subproblem data
        set_up_data(simdir)
    end
end

# This is the wrapper function to run all of parallelized benders
function parallelized_ptdf_benders(superdir)

    # Read the full config file
    config_file = joinpath(superdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    dates = toml_data["dates"]
    rep_prob = toml_data["representative_prob"]

    date_weights = Dict{Int, Tuple{String, Float64}}()
    for i in 1:length(dates)
        date_weights[i] = (dates[i], rep_prob[i])
    end

    # Set up master and subproblem directories
    set_up_directories(superdir)

    # Set up the data by reading the config
    master_data = set_up_data(superdir)

    # Compute the initial core point
    compute_superset_core_point(superdir)
    println("[DEBUG] Finished setting up directories and files for $superdir")

    # Create master problem
    master, y, theta = load_or_make_master_problem(superdir, master_data, date_weights)

    # Main benders loop
    gamma_val, s_power_val, s_energy_val = master_benders_loop(superdir, master, y, theta, master_data)

    # Save final investments
    export_investments_csv(master_data, gamma_val, s_power_val, s_energy_val, output_dir=joinpath(superdir,"final_output"))

end

# This function will get the initial core point for the benders (an over-invested first stage)
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

function master_benders_loop(superdir, master, y, theta, master_data, max_iterations=100000, tolerance=0.01)
    converged = false
    iter = master.ext[:iter]
    date_weights = master.ext[:date_weights]
    
    # Unpack y variables
    gamma, s_power, s_energy = y
    gamma_val, s_power_val, s_energy_val = nothing, nothing, nothing

    while !converged && iter < max_iterations
        # Solve master problem
        println("[DEBUG] Solving master problem: iteration $iter")
        update_master_objective!(superdir, master, master_data, y, theta)
        optimize!(master)
        
        if !has_values(master)
            error("Master problem has no solution: $(termination_status(master))")
        end
        
        # Get master solution
        gamma_val = value.(gamma)
        s_power_val = value.(s_power)
        s_energy_val = value.(s_energy)
        theta_val = value.(theta)
        y_raw = [gamma_val, s_power_val, s_energy_val]
        push!(master.ext[:y_raw], y_raw)

        # Save raw solution
        export_investments_csv(master_data, gamma_val, s_power_val, s_energy_val, output_dir=joinpath(superdir,"benders_output"), file_suffix="$iter")

        # Compute the evaluation point where Bender's cuts will be made
        y_eval, y_core = compute_eval_core_points(superdir, master, master_data)
        lambda_val = master.ext[:stabilization_lambda][end]

        # This is thread SAFE - pre-allocate and use indexed assignment
        date_weights_vec = sort(collect(date_weights), by=first) 
        subproblem_results = Vector{Tuple}(undef, length(date_weights_vec))
        @threads for i in 1:length(date_weights_vec)
            rep_index, (a_date, prob) = date_weights_vec[i]
            @assert rep_index == i  # Verify our assumption
            thread_id = threadid()
            println("Thread $thread_id: Starting subproblem for $a_date (rep_index $rep_index)")
            start_time = time()
            
            try
                simdir = joinpath(superdir, a_date)
                _, phi_val, duals = solve_subproblem_ptdf(simdir, y_eval, logging="Eval Point iter $iter for $a_date")
                
                solve_time = time() - start_time
                println("Thread $thread_id: Completed $a_date in $(round(solve_time, digits=2)) seconds")
                
                subproblem_results[i] = (a_date, simdir, phi_val, duals)
                
            catch e
                solve_time = time() - start_time
                println("Thread $thread_id: ERROR solving $a_date after $(round(solve_time, digits=2)) seconds: $e")
                rethrow(e)
            end
        end
        println("[DEBUG] Iteration $(iter): all subproblems solved!")

        # Sequentially save subproblem Bender's logs and add cuts to the master problem
        master_obj_val = objective_value(master)
        for rep_index in 1:length(subproblem_results)
            a_date, simdir, phi_val, duals = subproblem_results[rep_index]
            
            # Individual subproblem Bender's log
            benders_ptdf_write_to_csv(simdir, master_obj_val, theta_val[rep_index], phi_val, y_eval, lambda_val=lambda_val)

            # Add the cut to the master problem
            add_benders_cut_ptdf(master, theta[rep_index], duals, y, y_eval, phi_val)
        end

        # Master problem Bender's log
        R = length(date_weights)
        total_theta_val = sum(theta_val[r] * date_weights[r][2] for r in 1:R)
        total_phi_val = sum(subproblem_results[r][3] * date_weights[r][2] for r in 1:R)
        benders_ptdf_write_to_csv(superdir, master_obj_val, total_theta_val, total_phi_val, y_eval, lambda_val=lambda_val)

        # Clear references
        subproblem_results = nothing
        GC.gc()

        # Check convergence
        gap = abs(total_theta_val - total_phi_val) / (1e-10 + abs(total_theta_val))
        push!(master.ext[:gap], gap)
        println("[DEBUG] Benders iteration $(iter): Master objective = $(master_obj_val), theta = $(total_theta_val), phi = $(total_phi_val), gap = $(gap)")

        # If converged, exit benders loop
        lambda = master.ext[:stabilization_lambda][end]
        if gap < 0.01 && lambda < 1e-10
            println("[DEBUG] Final gap at discrete solution = $(gap) < $(tolerance). Converged after $(iter) iterations.")
            println("[DEBUG] Stabilization lambda was $(lambda).")
            converged = true
        end
        
        # Update iter count and save master problem with metadata/extensions
        master.ext[:iter] += 1
        iter = master.ext[:iter]
        save_master_problem(master, superdir)
    end

    if !converged
        println("[DEBUG] Benders decomposition did not converge after $max_iterations iterations")
    end
    
    # Return final solution
    return [gamma_val, s_power_val, s_energy_val]
end

    
    

