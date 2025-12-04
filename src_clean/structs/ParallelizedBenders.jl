using JuMP
using Gurobi
using CSV, DataFrames
using Base.Threads
using JSON
using TOML

include("../helpers/helpers.jl")
include("../utils/utils.jl")
include("BendersSubproblem.jl")
include("BendersMasterProblem.jl")

# Persistent subproblem worker for parallel execution
mutable struct PersistentSubproblemWorker
    subproblem::BendersSubproblem
    work_channel::Channel{Vector{Any}}
    result_channel::Channel{Any}
    task::Task
    worker_id::Int
    date::String
    simdir::String
end

# Main parallelized Benders decomposition struct
mutable struct ParallelizedBenders
    master::BendersMasterProblem
    workers::Vector{PersistentSubproblemWorker}
    superdir::String
    date_weights::Dict{Int, Tuple{String, Float64}}
    
    # Convergence parameters
    max_iterations::Int
    tolerance::Float64
    converged::Bool
    
    # Statistics tracking
    master_times::Vector{Float64}
    iteration_times::Vector{Float64}
    gaps::Vector{Float64}
    objectives::Vector{Float64}
    
    function ParallelizedBenders(superdir::String; 
                                max_iterations::Int=100000,
                                tolerance::Float64=0.005)
        
        # Read configuration
        config_file = joinpath(superdir, "config.toml")
        toml_data = TOML.parsefile(config_file)
        dates = toml_data["dates"]
        rep_prob = toml_data["representative_prob"]
        
        # Create date_weights mapping
        date_weights = Dict{Int, Tuple{String, Float64}}()
        for i in 1:length(dates)
            date_weights[i] = (dates[i], rep_prob[i])
        end
        
        # Set up directories and data
        set_up_directories(superdir, toml_data)
        master_data = set_up_data(superdir)
        
        println("[DEBUG] Finished setting up directories and files for $superdir")
        
        # Create master problem
        master = BendersMasterProblem(superdir, master_data, date_weights)
        
        # Initialize workers
        workers = create_persistent_workers(superdir, date_weights)
        
        # Let workers initialize
        sleep(1)
        
        new(master, 
            workers, 
            superdir, 
            date_weights,
            max_iterations, 
            tolerance, 
            false, # converged
            Float64[], # master_times
            Float64[], # iteration_times
            Float64[], # gaps
            Float64[] # objectives
            )
    end
end

# Set up directories for master and subproblems
function set_up_directories(superdir::String, toml_data::Dict)
    """
    Set up the superdirectory and subproblem directories, and generate data.json files.
    """
    dates = toml_data["dates"]
    
    # Create master output directories
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
        
        # Create subproblem data
        set_up_data(simdir)
    end
end

# Create persistent worker threads
function create_persistent_workers(superdir::String, 
                                  date_weights::Dict{Int, Tuple{String, Float64}})
    """
    Create persistent worker threads for parallel subproblem solving.
    Each worker maintains its own BendersSubproblem (wrapping PTDFModel).
    """
    workers = PersistentSubproblemWorker[]
    
    for i in 1:length(date_weights)
        date = date_weights[i][1]
        simdir = joinpath(superdir, date)
        work_channel = Channel{Vector{Any}}(1)
        result_channel = Channel{Any}(1)
        subproblem_data = JSON.parsefile(joinpath(simdir, "data.json"))
        
        # Create placeholder - will be initialized in worker thread
        worker = PersistentSubproblemWorker(
            BendersSubproblem(subproblem_data, Gurobi.Optimizer, simdir),  # Placeholder
            work_channel,
            result_channel,
            Task(() -> nothing),
            i,
            date,
            simdir
        )
        
        worker.task = Threads.@spawn worker_loop(worker)
        
        push!(workers, worker)
    end
    
    return workers
end

# Main worker loop running in separate thread
function worker_loop(worker::PersistentSubproblemWorker)
    """
    Main loop for a persistent worker thread.
    Initializes BendersSubproblem once, then repeatedly solves with different investments.
    """
    try
        # Load data
        data = JSON.parsefile(joinpath(worker.simdir, "data.json"))
        
        # Create BendersSubproblem (which creates PTDFModel internally)
        worker.subproblem = BendersSubproblem(
            data, 
            Gurobi.Optimizer, 
            worker.simdir;
            max_ptdf_iterations=256,
            max_ptdf_per_iteration=32,
            ptdf_tol=1e-6
        )
        
        println("[DEBUG] Worker $(worker.worker_id) ($(worker.date)) initialized on thread $(Threads.threadid())")
        
        # Main work loop
        while true
            # Wait for work (y_val)
            y_val = take!(worker.work_channel)
            
            # Solve subproblem
            result = solve_worker_subproblem(worker, y_val, data)
            
            # Return result
            put!(worker.result_channel, result)
        end
        
    catch e
        if e isa InvalidStateException
            println("[DEBUG] Worker $(worker.worker_id) shutting down gracefully")
        else
            println("[ERROR] Worker $(worker.worker_id) encountered error: $e")
            rethrow(e)
        end
    end
end

# Solve subproblem for a worker
function solve_worker_subproblem(worker::PersistentSubproblemWorker, 
                                y_val::Vector, 
                                data::Dict)
    """
    Solve the subproblem for given investment decisions.
    Returns a dictionary with results including duals for cut generation.
    """
    gamma_val, s_energy_val = y_val
    
    # Fix investment decisions
    fix_investments!(worker.subproblem, gamma_val, s_energy_val)
    
    # Phase 1: Check feasibility (minimize unserved energy)
    solve!(worker.subproblem, :feasibility; configure_optimizer=false)
    total_ue = get_total_unserved_energy(worker.subproblem)
    duals = extract_duals(worker.subproblem)
    
    if total_ue > 1e-6
        # Infeasible - need feasibility cut
        println("[DEBUG] Worker $(worker.worker_id) ($(worker.date)): FEASIBILITY cut, unserved energy = $(round(total_ue, digits=4))")
        
        return Dict(
            :total_ue => total_ue,
            :phi_val => 0.0,
            :duals => duals,
            :worker_id => worker.worker_id,
            :is_feasibility_cut => true,
            :date => worker.date
        )
    else
        # Feasible - optimize operations
        solve!(worker.subproblem, :operational; configure_optimizer=false)
        phi_val = get_objective_value(worker.subproblem)
        duals = extract_duals(worker.subproblem)
        total_ue_actual = get_total_unserved_energy(worker.subproblem)
        
        println("[DEBUG] Worker $(worker.worker_id) ($(worker.date)): OPTIMALITY cut, phi = $(round(phi_val, digits=2))")
        
        return Dict(
            :total_ue => total_ue_actual,
            :phi_val => phi_val,
            :duals => duals,
            :worker_id => worker.worker_id,
            :is_feasibility_cut => false,
            :date => worker.date
        )
    end
end

# Dispatch work to all workers in parallel
function solve_subproblems_parallel!(benders::ParallelizedBenders, y_val::Tuple)
    """
    Dispatch investment decisions to all workers and collect results in parallel.
    """
    gamma_val, s_energy_int_val = y_val
    y_vec = [gamma_val, s_energy_int_val]
    
    # Send work to all workers (concurrent dispatch)
    for worker in benders.workers
        put!(worker.work_channel, copy(y_vec))
    end
    
    # Collect all results (blocks until all workers are done)
    all_results = []
    for worker in benders.workers
        result = take!(worker.result_channel)
        push!(all_results, result)
    end
    
    return all_results
end

# Add all Benders cuts from parallel results
function add_all_cuts!(benders::ParallelizedBenders, 
                      theta_val::Vector{Float64},
                      all_results::Vector,
                      y_val::Tuple)
    """
    Add all Benders cuts from parallel subproblem solves.
    Returns aggregate statistics for convergence checking.
    """
    master = benders.master
    master_obj = get_objective_value(master)
    
    total_theta_val = 0.0
    total_phi_val = 0.0
    total_ue = 0.0
    
    for result in all_results
        ue = result[:total_ue]
        phi_val = result[:phi_val]
        duals = result[:duals]
        worker_id = result[:worker_id]
        is_feasibility_cut = result[:is_feasibility_cut]
        date = result[:date]
        
        if is_feasibility_cut
            # Add feasibility cut
            add_feasibility_cut!(master, duals, y_val, ue)
        else
            # Add optimality cut for this representative period
            add_benders_cut!(master, worker_id, duals, y_val, phi_val)
        end
        
        # Create rep-day specific log
        output_dir = joinpath(benders.superdir, date)
        benders_ptdf_write_to_csv(output_dir, y_val, master_obj, 
                                 theta_val[worker_id], phi_val, ue)
        
        # Accumulate for aggregate statistics
        weight = benders.date_weights[worker_id][2]
        total_theta_val += theta_val[worker_id] * weight
        
        if is_feasibility_cut
            penalty = master.data["param"]["under_served_penalty"]
            total_phi_val += (theta_val[worker_id] + penalty * ue) * weight
        else
            total_phi_val += phi_val * weight
        end
        
        total_ue += ue * weight
    end
    
    # Write aggregate log
    benders_ptdf_write_to_csv(benders.superdir, y_val, master_obj, 
                             total_theta_val, total_phi_val, total_ue)
    
    current_obj = master_obj + total_phi_val - total_theta_val
    
    return current_obj, total_phi_val, total_theta_val, total_ue
end

# Update trust region based on iteration results
function update_trust_region!(benders::ParallelizedBenders,
                             y_val::Tuple,
                             current_obj::Float64,
                             total_theta_val::Float64,
                             total_phi_val::Float64,
                             total_ue::Float64)
    """
    Update trust region constraints based on iteration results.
    Implements serious step vs null step logic.
    """
    master = benders.master
    
    if master.stabilization != "trust_region"
        return
    end
    
    if total_ue < 1e-6
        if current_obj < master.upper_bound
            # Actual improvement - serious step
            println("[DEBUG] Serious step: objective improved from $(master.total_obj[end]) to $current_obj")
            push!(master.y_trust, [y_val[1], y_val[2]])
            master.upper_bound = current_obj
            
            # Reset l1 radius
            println("[DEBUG] Resetting l1_radius to 1")
            master.jump_model.ext[:l1_radius] = 
                vcat(get(master.jump_model.ext, :l1_radius, Int[]), [1])
            
            # Add level set
            add_level_set!(master, current_obj)
        else
            # No improvement - null step
            l1_distance = compare_y_vals([y_val[1], y_val[2]], master.last_y_val)
            
            if isapprox(total_theta_val, total_phi_val, atol=1e-3) || l1_distance == 0
                # Stuck in local region - expand trust region
                current_radius = get(master.jump_model.ext, :l1_radius, [0])[end]
                new_radius = current_radius + 1
                println("[DEBUG] Null step: expanding l1_radius from $current_radius to $new_radius")
                master.jump_model.ext[:l1_radius] = 
                    vcat(get(master.jump_model.ext, :l1_radius, Int[]), [new_radius])
            end
        end
    end
end

# Check convergence
function check_convergence(benders::ParallelizedBenders,
                          gap::Float64,
                          total_ue::Float64)
    """
    Check if Benders decomposition has converged.
    """
    if gap < benders.tolerance && total_ue < 1e-6
        benders.converged = true
        return true
    end
    return false
end

# Main solve function
function solve!(benders::ParallelizedBenders)
    """
    Execute the Benders decomposition algorithm.
    """
    master = benders.master
    
    try
        for iteration in 1:benders.max_iterations
            iter_start = time()
            
            println("\n" * "="^80)
            println("[ITERATION $(master.iter)]")
            println("="^80)
            
            # Solve master problem
            println("[DEBUG] Solving master problem...")
            add_trust_region!(master)
            master_start = time()
            solve!(master)
            master_time = time() - master_start
            push!(benders.master_times, master_time)
            println("[DEBUG] Master solve time: $(round(master_time, digits=2))s")
            
            # Get master solution
            gamma_val, s_energy_int_val = get_investments(master)
            theta_val = get_theta_values(master)
            master_obj = get_objective_value(master)
            y_val = (gamma_val, s_energy_int_val)
            
            println("[DEBUG] Master objective: $(round(master_obj, digits=2))")
            
            # Export investments
            export_investments_csv(master.data, gamma_val, s_energy_int_val,
                                 output_dir=joinpath(benders.superdir, "benders_output"),
                                 file_suffix="$(master.iter)")
            
            # Solve all subproblems in parallel
            println("[DEBUG] Solving $(length(benders.workers)) subproblems in parallel...")
            all_results = solve_subproblems_parallel!(benders, y_val)
            
            # Add Benders cuts
            println("[DEBUG] Adding Benders cuts...")
            current_obj, total_phi_val, total_theta_val, total_ue = 
                add_all_cuts!(benders, theta_val, all_results, y_val)
            
            # Calculate gap
            gap = calculate_gap(benders, master_obj, total_phi_val, total_theta_val)
            push!(benders.gaps, gap)
            push!(benders.objectives, current_obj)
            
            println("[DEBUG] Current objective: $(round(current_obj, digits=2))")
            println("[DEBUG] Gap: $(round(gap * 100, digits=4))%")
            println("[DEBUG] Total unserved energy: $(round(total_ue, digits=6))")
            
            # Update trust region
            update_trust_region!(benders, y_val, current_obj, 
                               total_theta_val, total_phi_val, total_ue)
            
            # Update tracking
            update_tracking!(master, y_val, total_ue, current_obj)
            
            # Record iteration time
            iter_time = time() - iter_start
            push!(benders.iteration_times, iter_time)
            println("[DEBUG] Iteration time: $(round(iter_time, digits=2))s")
            
            # Check convergence
            if check_convergence(benders, gap, total_ue)
                println("\n" * "="^80)
                println("[SUCCESS] Benders converged!")
                println("Final gap: $(round(gap * 100, digits=4))%")
                println("Final objective: $(round(current_obj, digits=2))")
                println("Total iterations: $(master.iter)")
                println("Total time: $(round(sum(benders.iteration_times), digits=2))s")
                println("="^80)
                export_results(benders)
                break
            end

            # Export results
            export_results(benders)
            
            # Increment iteration
            increment_iteration!(master)
        end
        
        if !benders.converged
            println("\n" * "="^80)
            println("[WARNING] Maximum iterations reached without convergence")
            println("Final gap: $(round(benders.gaps[end] * 100, digits=4))%")
            println("="^80)
        end
        
    finally
        # Clean shutdown of workers
        shutdown!(benders)
    end
    
    # Return final solution
    return master.y_trust[end]
end

function calculate_gap(benders::ParallelizedBenders, master_obj, total_phi_val, total_theta_val)
    """
    Calculates the gap of Benders: (UB - LB) / UB
    UB should be the lowest observed UB.
    LB should be the highest observed LB.
    """
    if benders.master.stabilization == "trust_region"
        gap = abs(benders.master.upper_bound - benders.master.lower_bound) /
                    (benders.master.upper_bound)

    else
        # vanilla calculation of benders gap (master_obj )
        gap = abs(benders.master.upper_bound - master_obj) / 
                  (benders.master.upper_bound)
    end

    return gap
end

# Shutdown workers
function shutdown!(benders::ParallelizedBenders)
    """
    Gracefully shutdown all worker threads.
    """
    println("\n[DEBUG] Shutting down $(length(benders.workers)) workers...")
    
    for worker in benders.workers
        try
            close(worker.work_channel)
        catch e
            println("[WARNING] Error closing channel for worker $(worker.worker_id): $e")
        end
    end
    
    # Wait for graceful shutdown
    sleep(1)
    println("[DEBUG] Shutdown complete")
end

# Export final results
function export_results(benders::ParallelizedBenders)
    """
    Export final investment decisions and statistics.
    """
    gamma_val, s_energy_int_val = benders.master.y_trust[end]
    
    # Export investments
    export_investments_csv(benders.master.data, gamma_val, s_energy_int_val,
                          output_dir=joinpath(benders.superdir, "output"))
    
    # Export convergence statistics
    stats_df = DataFrame(
        iteration = 1:length(benders.gaps),
        gap = benders.gaps,
        objective = benders.objectives,
        time = benders.iteration_times,
        master_time = benders.master_times
    )
    
    CSV.write(joinpath(benders.superdir, "output", "convergence_stats.csv"), stats_df)
    
    println("[DEBUG] Results exported to $(benders.superdir)")
end

# Main entry point
function parallelized_ptdf_benders(superdir::String; 
                                  max_iterations::Int=100000,
                                  tolerance::Float64=0.005)
    """
    Main entry point for parallelized Benders decomposition.
    """
    # TODO: make a wrapper function in main.jl
    # println("[DEBUG] Computing global lower bound...")
    # lb_value = compute_benders_lb(superdir, force_lb=force_lb)

    # Create Benders decomposition
    benders = ParallelizedBenders(superdir; 
                                 max_iterations=max_iterations,
                                 tolerance=tolerance)
    
    # Solve
    solve!(benders)
        
    return benders
end