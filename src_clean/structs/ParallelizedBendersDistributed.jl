using JuMP
using Gurobi
using CSV, DataFrames
using Distributed
using ClusterManagers
using JSON
using TOML

include("../helpers/helpers.jl")
include("../utils/utils.jl")
include("BendersSubproblem.jl")
include("BendersMasterProblem.jl")
include("MultistageBendersMasterProblem.jl")

# Worker initialization function - runs on each worker process
@everywhere function init_sub_worker(simdir::String, date::String)
    """
    Initialize a worker process with its subproblem.
    This function runs on each worker process.
    """
    # Load data
    data = set_up_data(simdir)
    
    # Create BendersSubproblem (which creates PTDFModel internally)
    subproblem = BendersSubproblem(
        data, 
        Gurobi.Optimizer, 
        simdir;
        max_ptdf_iterations=256,
        max_ptdf_per_iteration=32,
        ptdf_tol=1e-6
    )
    
    # Store in global scope on THIS worker (never serialize this!)
    global PERSISTENT_SUBPROBLEM = subproblem
    global PERSISTENT_DATA = data
    
    println("[DEBUG] Worker for $date initialized on process $(myid())")
    
    # Return ONLY non-Gurobi metadata
    return (date=date, simdir=simdir, initialized=true)
end

# Worker subproblem solve function - runs on each worker
@everywhere function solve_worker_subproblem(worker_state, y_val::Vector)
    """
    Solve using the persistent subproblem stored on this worker.
    Never serialize the Gurobi model - only return numerical results.
    """
    # Use the persistent subproblem (already on this worker)
    subproblem = PERSISTENT_SUBPROBLEM
    date = worker_state.date
    
    gamma_val, s_energy_val = y_val
    
    # Fix investment decisions
    fix_investments!(subproblem, gamma_val, s_energy_val)
    
    # Solve feasibility (PTDF cuts accumulate in persistent model)
    solve!(subproblem, :feasibility; configure_optimizer=false)
    total_ue = get_total_unserved_energy(subproblem)
    duals = extract_duals(subproblem)
    
    if total_ue > 1e-6
        println("[DEBUG] Worker $date (pid=$(myid())): FEASIBILITY cut, UE = $(round(total_ue, digits=4))")
        
        return Dict(
            :total_ue => total_ue,
            :phi_val => 0.0,
            :duals => duals,
            :is_feasibility_cut => true,
            :date => date
        )
    else
        solve!(subproblem, :operational; configure_optimizer=false)
        phi_val = get_objective_value(subproblem)
        duals = extract_duals(subproblem)
        total_ue_actual = get_total_unserved_energy(subproblem)
        
        println("[DEBUG] Worker $date (pid=$(myid())): OPTIMALITY cut, phi = $(round(phi_val, digits=2))")
        
        return Dict(
            :total_ue => total_ue_actual,
            :phi_val => phi_val,
            :duals => duals,
            :is_feasibility_cut => false,
            :date => date
        )
    end
end

# Main parallelized Benders decomposition struct
mutable struct ParallelizedBendersDistributed
    master::MultistageBendersMasterProblem
    worker_pids::Vector{Int}
    worker_dates::Vector{String}
    worker_states::Dict{String, Any}  # Store worker states by date
    superdir::String
    date_weights::Dict{String, Float64}
    
    years::Vector{Int}
    
    # Convergence parameters
    max_iterations::Int
    tolerance::Float64
    converged::Bool
    
    # Statistics tracking
    master_times::Vector{Float64}
    iteration_times::Vector{Float64}
    gaps::Vector{Float64}
    objectives::Vector{Float64}
    
    function ParallelizedBendersDistributed(superdir::String; 
                                max_iterations::Int=100000,
                                tolerance::Float64=0.005,
                                lb::Float64=0.0,
                                use_slurm::Bool=true)
        
        # Read configuration
        config_file = joinpath(superdir, "config.toml")
        toml_data = TOML.parsefile(config_file)
        dates = toml_data["dates"]
        rep_prob = toml_data["representative_prob"]
        years = get(toml_data, "years", [toml_data["decarbonization_year"]])

        # Create date_weights mapping
        date_weights = Dict{String, Float64}()
        for i in 1:length(dates)
            date_weights[dates[i][6:end]] = rep_prob[i]
        end

        # Clear directories
        clear_directory(joinpath(superdir, "benders_output"))
        clear_directory(joinpath(superdir, "output"))

        # Set up directories and data
        simdirs = create_simdirs(superdir, toml_data)
        master_data = set_up_data(superdir)
        
        println("[DEBUG] Finished setting up directories and files for $superdir")
        
        # Initialize master problem
        master = MultistageBendersMasterProblem(superdir, master_data, date_weights)
        
        master.lower_bound = lb
        
        # Add worker processes and initialize subproblems
        worker_pids, worker_dates, worker_states = setup_workers(superdir, simdirs, 
                                                                  length(simdirs), use_slurm)
        
        println("[DEBUG] Successfully initialized $(length(worker_pids)) workers")
        
        new(master, 
            worker_pids,
            worker_dates,
            worker_states,  # NEW: store worker states
            superdir, 
            date_weights,
            years,
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

# Setup worker processes
function setup_workers(superdir::String, simdirs, n_workers::Int, use_slurm::Bool)
    """
    Setup worker processes using either SLURM or local processes.
    One worker per representative day.
    Returns arrays of worker PIDs, dates, and worker states.
    """
    # n_workers should equal length(simdirs)
    n_workers = length(simdirs)
    
    println("[DEBUG] Setting up $n_workers workers (one per representative day)...")
    
    if use_slurm
        # Get SLURM environment variables
        slurm_ntasks = get(ENV, "SLURM_NTASKS", nothing)
        slurm_jobid = get(ENV, "SLURM_JOB_ID", nothing)
        
        if slurm_ntasks !== nothing
            println("[DEBUG] Running on SLURM with $slurm_ntasks tasks (Job ID: $slurm_jobid)")
            # Use SlurmManager to add workers
            addprocs(SlurmManager(n_workers))
        else
            @warn "use_slurm=true but SLURM environment not detected. Using local processes."
            addprocs(n_workers)
        end
    else
        # Use local processes
        println("[DEBUG] Adding $n_workers local worker processes...")
        addprocs(n_workers)
    end
    
    # Get worker PIDs
    worker_pids = workers()
    println("[DEBUG] Worker PIDs: $worker_pids")
    
    # Get absolute path to project directory
    project_dir = pwd()  # Should be /storage/home/hcoda1/1/kwu381/TNEP-Storage
    
    # Load code and packages on all workers
    println("[DEBUG] Loading packages and code on all workers...")
    load_worker_code(project_dir)
    println("[DEBUG] Code loaded on all workers")
    
    # Initialize one worker per simdir IN PARALLEL
    worker_dates = String[]
    worker_states = Dict{String, Any}()

    futures = Future[]
    init_data = []

    for (i, simdir) in enumerate(simdirs)
        date = basename(simdir)
        push!(worker_dates, date)
        simdir_full = joinpath(superdir, date)
        
        # Assign worker i to simdir i (one-to-one mapping)
        worker_pid = worker_pids[i]
        
        # Queue initialization (don't wait yet)
        println("[DEBUG] Queuing initialization for worker $(worker_pid) for date $date")
        future = remotecall(init_sub_worker, worker_pid, simdir_full, date)
        push!(futures, future)
        push!(init_data, (worker_pid, date))
    end

    # Wait for all initializations to complete
    println("[DEBUG] Waiting for all workers to initialize...")
    results = fetch.(futures)

    # Store results
    for (i, result) in enumerate(results)
        worker_pid, date = init_data[i]
        worker_states[date] = result
        println("[DEBUG] Worker $worker_pid initialized for date: $date")
    end

    return worker_pids, worker_dates, worker_states
end

# Helper function to load code on workers (runtime version)
function load_worker_code(project_dir::String)
    """Load packages and source files on all worker processes in parallel."""
    @sync begin
        for w in workers()
            @async begin
                result = remotecall_fetch(w) do
                    Core.eval(Main, :(using Pkg))
                    Core.eval(Main, :(Pkg.activate($project_dir)))
                    Core.eval(Main, :(cd($project_dir)))
                    Core.eval(Main, :(include(joinpath($project_dir, "src_clean/main.jl"))))
                    return "Worker $(myid()) loaded successfully"
                end
                println("[DEBUG] $result")
            end
        end
    end
end

# Dispatch work to all workers in parallel
function solve_subproblems_parallel!(benders::ParallelizedBendersDistributed, y_val)
    """
    Dispatch investment decisions to all workers and collect results in parallel.
    One worker per date (one-to-one mapping).
    
    Multistage:   y_val is Dict{Int => Vector{Vector, Vector}}
                  Send year-specific investments to each worker based on date
    """
    # Create futures for all workers (one-to-one mapping)
    futures = Future[]
    
    for (i, date) in enumerate(benders.worker_dates)
        # Get the worker state for this date
        worker_state = benders.worker_states[date]
        
        # Use direct one-to-one mapping: worker i handles date i
        worker_pid = benders.worker_pids[i]
        
        # Get year and prepare y_vec
        # Multi-stage: send year-specific investments
        year = parse(Int, date[begin:4])
        y_vec = collect(y_val[year])
        
        # Call solve_worker_subproblem on the remote worker
        future = remotecall(solve_worker_subproblem, worker_pid, 
                           worker_state, y_vec)
        push!(futures, future)
    end
    
    # Collect all results (blocks until all workers are done)
    all_results = fetch.(futures)
    
    return all_results
end

# Add all Benders cuts from parallel results
function add_all_cuts!(benders::ParallelizedBendersDistributed, 
                      theta_val,
                      all_results::Vector,
                      y_val)
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
        is_feasibility_cut = result[:is_feasibility_cut]
        date = result[:date]
        year = parse(Int, date[begin:4])
        
        if is_feasibility_cut
            add_feasibility_cut!(master, year, duals, y_val[year], ue)
        else
            add_benders_cut!(master, date, year, duals, y_val[year], phi_val)
        end
        
        # Create rep-day specific log
        output_dir = joinpath(benders.superdir, date)
        benders_ptdf_write_to_csv(output_dir, y_val[year], master_obj, 
                                    theta_val[date], phi_val, ue)
        
        # Accumulate for aggregate statistics
        weight = benders.date_weights[date[6:end]]
        disc_factors = master.disc_factors
        
        total_theta_val += disc_factors[year] * theta_val[date] * weight
        
        if is_feasibility_cut
            penalty = master.data["param"]["under_served_penalty"]
            total_phi_val += disc_factors[year] * (theta_val[date] + penalty * ue) * weight
        else
            total_phi_val += disc_factors[year] * phi_val * weight
        end
        
        total_ue += ue * weight
    end
    
    benders_ptdf_write_to_csv(benders.superdir, y_val[master.years[1]], master_obj, 
                            total_theta_val, total_phi_val, total_ue)

    current_obj = master_obj + total_phi_val - total_theta_val
    
    return current_obj, total_phi_val, total_theta_val, total_ue
end

# Update trust region based on iteration results
function update_trust_region!(benders::ParallelizedBendersDistributed,
                             y_val,
                             current_obj::Float64,
                             total_theta_val::Float64,
                             total_phi_val::Float64,
                             total_ue::Float64)
    """
    Update trust region constraints based on iteration results.
    """
    master = benders.master
    
    if master.stabilization != "trust_region"
        return
    end
    
    if total_ue < 1e-6
        if current_obj < master.upper_bound
            println("[DEBUG] Serious step: objective improved from $(master.total_obj[end]) to $current_obj")
            push!(master.y_trust, y_val)
            master.upper_bound = current_obj

            println("[DEBUG] Resetting l1_radius to 1")
            master.jump_model.ext[:l1_radius] = 
                vcat(get(master.jump_model.ext, :l1_radius, Int[]), [1])
            
            add_level_set!(master, current_obj)
        else
            println("[DEBUG] Null step")
            if all(compare_y_vals(y_val[k], master.last_y_val[k]) == 0 for k in keys(y_val))
                current_radius = get(master.jump_model.ext, :l1_radius, [0])[end]
                new_radius = current_radius + 1
                println("[DEBUG] Repeat solution: expanding l1_radius from $current_radius to $new_radius")
                master.jump_model.ext[:l1_radius] = 
                    vcat(get(master.jump_model.ext, :l1_radius, Int[]), [new_radius])
            end
        end
    end
end

# Check convergence
function check_convergence(benders::ParallelizedBendersDistributed,
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
function solve!(benders::ParallelizedBendersDistributed)
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
            result = solve!(master)
            if result !== nothing
                @error "Master infeasible - returning for debugging"
                return (converged=false, master=master)
            end
            master_time = time() - master_start
            push!(benders.master_times, master_time)
            println("[DEBUG] Master solve time: $(round(master_time, digits=2))s")
            
            y_val = get_investments(master)
            for year in master.years
                export_investments_csv(master.data, y_val[year][1], y_val[year][2],
                                    output_dir=joinpath(benders.superdir, "benders_output"),
                                    file_suffix="$(master.iter)_$(year)")
            end
            
            theta_val = get_theta_values(master)
            master_obj = get_objective_value(master)
            
            println("[DEBUG] Master objective: $(round(master_obj, digits=2))")
            
            # Solve all subproblems in parallel
            println("[DEBUG] Solving $(length(benders.worker_dates)) subproblems across $(length(benders.worker_pids)) workers...")
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

    return master.y_trust[end]
end

function calculate_gap(benders::ParallelizedBendersDistributed, master_obj, total_phi_val, total_theta_val)
    """
    Calculates the gap of Benders: (UB - LB) / UB
    """
    if benders.master.stabilization == "trust_region"
        gap = abs(benders.master.upper_bound - benders.master.lower_bound) /
                    (benders.master.upper_bound)
    else
        gap = abs(benders.master.upper_bound - master_obj) / 
                  (benders.master.upper_bound)
    end
    return gap
end

# Shutdown workers
function shutdown!(benders::ParallelizedBendersDistributed)
    """
    Gracefully shutdown all worker processes.
    """
    println("\n[DEBUG] Shutting down $(length(benders.worker_pids)) worker processes...")
    
    # Remove worker processes
    rmprocs(benders.worker_pids)
    
    println("[DEBUG] Shutdown complete")
end

# Export final results
function export_results(benders::ParallelizedBendersDistributed)
    """
    Export final investment decisions and statistics.
    """
    y_val = benders.master.y_trust[end]
    for year in benders.master.years
        ans = y_val[year]
        export_investments_csv(benders.master.data, ans[1], ans[2],
                        output_dir=joinpath(benders.superdir, "output"), file_suffix="$(year)")
    end
    
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
function parallelized_ptdf_benders_distributed(superdir::String; 
                                  max_iterations::Int=100000,
                                  tolerance::Float64=0.005,
                                  lb::Float64=0.0,
                                  use_slurm::Bool=true)
    """
    Main entry point for parallelized Benders decomposition using Distributed.
    
    Arguments:
    - superdir: Directory containing problem data
    - max_iterations: Maximum number of Benders iterations
    - tolerance: Convergence tolerance
    - lb: Lower bound estimate
    - use_slurm: Whether to use SLURM cluster manager
    """
    # Create Benders decomposition
    benders = ParallelizedBendersDistributed(superdir; 
                                 max_iterations=max_iterations,
                                 tolerance=tolerance,
                                 lb=lb,
                                 use_slurm=use_slurm)
    
    # Solve
    result = solve!(benders)

    # Check if we returned early due to infeasibility
    if result isa NamedTuple && haskey(result, :converged) && !result.converged
        @warn "Returning early due to master infeasibility"
        result[2] = master
        model = master.jump_model
        print_conflict!(model)        
    end
        
    return benders
end

"""
    print_conflict!(model)
Compute and print a conflict for an infeasible `model`.
"""
function print_conflict!(model)
    JuMP.compute_conflict!(model)
    ctypes = list_of_constraint_types(model)
    for (F, S) in ctypes
        cons = all_constraints(model, F, S)
        for i in eachindex(cons)
            isassigned(cons, i) || continue
            con = cons[i]
            cst = MOI.get(model, MOI.ConstraintConflictStatus(), con)
            cst == MOI.IN_CONFLICT && @info name(con) con
        end
    end
    return nothing
end
