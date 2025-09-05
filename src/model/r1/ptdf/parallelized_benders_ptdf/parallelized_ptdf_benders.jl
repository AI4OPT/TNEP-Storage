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
    println("[DEBUG] Finished setting up directories and files for $superdir")

    # Create master problem
    master, y, theta = define_master_ptdf(superdir, master_data, date_weights)

    # Main benders loop
    gamma_val, s_energy_int_val = master_benders_loop(superdir, master, master_data, date_weights)

    # Save final investments
    export_investments_csv(master_data, gamma_val, s_energy_int_val, output_dir=joinpath(superdir,"final_output"))
    return master
end

mutable struct PersistentSubproblemWorker
    sub_model::Model
    work_channel::Channel{Vector{Any}} # receives y_val
    result_channel::Channel{Any} # returns duals
    task::Task
    worker_id::Int
end

function create_persistent_workers(superdir, date_weights)
    workers = PersistentSubproblemWorker[]

    for i in 1:length(date_weights)
        work_channel = Channel{Vector{Any}}(1)
        result_channel = Channel{Any}(1)

        worker = PersistentSubproblemWorker(
            Model(),
            work_channel,
            result_channel,
            Task(() -> nothing),
            i
        )

        worker.task = Threads.@spawn begin
            # Create subproblem model once in this thread
            id = worker.worker_id
            date = date_weights[id][1]
            simdir = joinpath(superdir, date)
            data = JSON.parsefile(joinpath(simdir, "data.json"))

            # Initialize tracked constraints dictionary to store across iterations
            tracked_constraints = Dict{Tuple{Int,Int,Int,Bool}, Bool}()

            # Load previous constraints if available
            if isfile(joinpath(simdir, "tracked_constraints.csv"))
                tracked_df = CSV.read(joinpath(simdir, "tracked_constraints.csv"), DataFrame)
                println("[DEBUG] Loading $(nrow(tracked_df)) previously tracked constraints")
                
                for row in eachrow(tracked_df)
                    if row.tracked
                        tracked_constraints[(row.arc, row.rep, row.time, row.ub)] = true
                    end
                end
            end

            worker.sub_model = create_subproblem_model(simdir, data, tracked_constraints)
            println("Worker $(worker.worker_id) on date $(date) initialized on thread $(Threads.threadid())")

            while true
                try
                    # Wait for work (y_val)
                    y_val = take!(work_channel)

                    # Solve subproblem
                    set_sub_investments!(worker.sub_model, y_val)
                    # Phase 1: check feasibility of no load shed
                    set_underserved_objective!(worker.sub_model)
                    solve_subproblem_ptdf!(worker.sub_model; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6)
                    total_ue = objective_value(worker.sub_model)
                    duals = extract_subproblem_duals(worker.sub_model)

                    if total_ue > 0 # If failed Phase 1, infeasible due to load shed
                        println("[DEBUG] Worker $(worker.worker_id) on date $(date): add FEASIBILITY cut, load shed $(total_ue) detected")

                        results = [total_ue, 0, duals, worker.worker_id]
                        put!(result_channel, results)
                    
                    else # Passed Phase 1, move onto Phase 2 to check operational objective
                        set_operational_objective!(worker.sub_model)
                        solve_subproblem_ptdf!(worker.sub_model; max_ptdf_iterations=256, max_ptdf_per_iteration=32, ptdf_tol=1e-6)
                        phi_val = objective_value(worker.sub_model)
                        duals = extract_subproblem_duals(worker.sub_model)
                        total_ue = vec(sum(value.(worker.sub_model[:ue]),dims=[1,2,3]))[1]

                        results = [total_ue, phi_val, duals, worker.worker_id]
                        put!(result_channel, results)
                    end
                
                catch e
                    if e isa InvalidStateException
                        println("Worker $(worker.worker_id) shutting down")
                        break  # Channel closed, exit gracefully
                    else
                        println("Error in worker $(worker.worker_id): $e")
                        put!(result_channel, nothing)  # Signal error
                    end
                end
            end
        end

        push!(workers, worker)
    end

    return workers
end

function solve_subproblems_parallel!(workers, y_val)
    # Send work to all workers (concurrent dispatch)
    for worker in workers
        put!(worker.work_channel, copy(y_val))
    end

    # Collect all results (wait for completion)
    all_results = []
    for worker in workers   
        results = take!(worker.result_channel) # Blocks until worker is done
        push!(all_results, results)
    end

    return all_results
end

function add_all_benders_cuts(master, master_obj, theta_val, all_results, y_val)
    y = master[:gamma], master[:s_energy_int]
    date_weights = master.ext[:date_weights]
    superdir = master.ext[:superdir]
    data = master.ext[:data]

    total_theta_val = 0
    total_phi_val = 0
    total_ue = 0

    for result in all_results
        ue = result[1]
        phi_val = result[2]
        duals = result[3]
        date_index = result[4]

        if ue > 0
            # Add feasibility cut
            add_benders_cut_ptdf(master, master[:ue_sum], duals, y, y_val, ue)
        else
            # Add optimality cut
            add_benders_cut_ptdf(master, master[:theta][date_index], duals, y, y_val, phi_val)
        end

        # Create the rep-day specific benders logs
        output_dir = joinpath(superdir, date_weights[date_index][1])
        benders_ptdf_write_to_csv(output_dir, y_val, master_obj, theta_val[date_index], phi_val, ue)

        # Also begin building the aggregate benders log
        total_theta_val += theta_val[date_index] * date_weights[date_index][2]
        if ue > 0
            total_phi_val += (theta_val[date_index] + data["param"]["under_served_penalty"] * ue) * date_weights[date_index][2]
        else
            total_phi_val += phi_val * date_weights[date_index][2]
        end
        total_ue += ue * date_weights[date_index][2]
    end

    benders_ptdf_write_to_csv(superdir, y_val, master_obj, total_theta_val, total_phi_val, total_ue)
    current_obj = master_obj + total_phi_val - total_theta_val
    return current_obj, total_ue
end

function master_benders_loop(superdir, master, master_data, date_weights, max_iterations=100000, tolerance=0.01)
    max_iterations=100000
    tolerance=0.01
    converged = false
    iter = master.ext[:iter]
    date_weights = master.ext[:date_weights]

    workers = create_persistent_workers(superdir, date_weights)
    # Let workers initialize
    sleep(1)

    gamma_val, s_energy_int_val = nothing, nothing

    try
        for iteration in 1:max_iterations
            println("[DEBUG] Solving master problem: iteration $iter")
            # Solve master problem
            add_trust_region!(master)
            optimize!(master)

            # Get master solution
            gamma_val = value.(master[:gamma])
            s_energy_int_val = value.(master[:s_energy_int])
            theta_val = value.(master[:theta])

            # Save lower bound
            master_obj = objective_value(master)

            # Export investments
            y_val = [gamma_val, s_energy_int_val]
            export_investments_csv(master_data, gamma_val, s_energy_int_val, output_dir=joinpath(superdir,"benders_output"), file_suffix="$iter")

            # Solve all subproblems in parallel
            println("[DEBUG] Iteration $(iter): solving all subproblems in parallel...")
            all_results = solve_subproblems_parallel!(workers, y_val)

            # Add Benders cuts sequentially, make Benders logs
            println("[DEBUG] Iteration $(iter): adding Benders cuts...")
            current_obj, total_ue = add_all_benders_cuts(master, master_obj, theta_val, all_results, y_val)

            # Update with serious step of the trust/anchor point if satisfy conditions
            if master.ext[:stabilization] == "trust_region"
                l1_max = get(master_data["param"], "l1_max", 8192)
                if total_ue > 0
                    new_l1_radius = minimum([l1_max, master.ext[:l1_radius][end] * 2])
                    println("[DEBUG] Load shed $total_ue detected, expanding transmission l1_radius to $(new_l1_radius)")
                    push!(master.ext[:l1_radius], new_l1_radius)
                else
                    if (current_obj < master.ext[:total_obj][end])
                        println("[DEBUG] Cheaper total objective $(current_obj) detected, taking serious step")
                        push!(master.ext[:y_trust], y_val)
                        push!(master.ext[:total_obj], current_obj)

                        new_l1_radius = maximum([1, master.ext[:l1_radius][end] / 2])
                        println("[DEBUG] Serious step, so reducing l1_radius to $(new_l1_radius)")
                        push!(master.ext[:l1_radius], new_l1_radius)
                    else
                        new_l1_radius = minimum([l1_max, master.ext[:l1_radius][end] * 2])
                        println("[DEBUG] No load shed or improvement in objective, expanding transmission l1_radius to $(new_l1_radius)")
                        push!(master.ext[:l1_radius], new_l1_radius)
                    end
                end
            end

            # TODO: check for convergence
            master.ext[:iter] += 1
            iter = master.ext[:iter]
        end

    finally
        # Clean shutdown
        println("Shutting down workers...")
        for worker in workers
            close(worker.work_channel)
        end

        # Wait a moment for graceful shutdown
        sleep(1)
    end

    # Return final solution
    return [gamma_val, s_energy_int_val]
end