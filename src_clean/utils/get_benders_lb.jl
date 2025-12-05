using TOML
using CSV
using DataFrames

include("../structs/ExpansionPlanner.jl")

function create_simdirs(superdir::String, toml_data::Dict; only_feasibility=nothing)
    simdirs = String[]
    year = toml_data["decarbonization_year"]
    
    for rep in toml_data["dates"]
        simdir = joinpath(superdir, string(year) * rep[5:end])
        
        # Create directory if needed
        mkpath(simdir)
        
        # Create modified config
        config = deepcopy(toml_data)
        config["dates"] = [rep]
        config["representative_prob"] = [1.0]
        config["num_representatives"] = 1

        if !isnothing(only_feasibility)
            config["only_feasibility"] = only_feasibility
        end
        
        # Write config
        open(joinpath(simdir, "config.toml"), "w") do io
            TOML.print(io, config)
        end
        
        push!(simdirs, simdir)
        println("Created simulation directory: $simdir")
    end
    
    return simdirs
end

function create_sbatch_file(job_name::String, 
                           log_file::String, 
                           sim_path::String, 
                           output_file::String,
                           stage_type::String;  # "first_stage" or "second_stage"
                           nodes::Int=1,
                           tasks_per_node::Int=24,
                           mem_per_cpu::String="12G",
                           time_limit::String="24:00:00",
                           partition::String="inferno",
                           account::String="gts-phentenryck3-coda20",
                           email::String="kwu381@gatech.edu")
    content = """
    #!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -q$partition
    #SBATCH --account=$account
    #SBATCH -N$nodes --ntasks-per-node=$tasks_per_node
    #SBATCH --mem-per-cpu=$mem_per_cpu
    #SBATCH -t$time_limit
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$log_file.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=$email
    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=. exp/clean/run_lb_single.jl $sim_path $tasks_per_node $stage_type
    """
    open(output_file, "w") do file
        write(file, content)
    end
    println("Created sbatch file: $output_file")
end

function create_batch_files(simdirs::Vector{String},
                           pace_dir::String,
                           name_suffix::String,
                           stage_type::String,
                           submit_jobs::Bool)
    mkpath(pace_dir)
    for simdir in simdirs
        rep = basename(simdir)
        job_name = "$(rep)$(name_suffix)"
        log_file = joinpath("lower_bound", job_name)
        output_file = joinpath(pace_dir, "run_$(rep).sbatch")
        create_sbatch_file(job_name, log_file, simdir, output_file, stage_type)
    end
    if submit_jobs
        submit_sbatch_jobs(pace_dir)
    end
end

function submit_sbatch_jobs(dir_name::String)
    run(`bash -c "for file in \"$dir_name\"/*.sbatch; do sbatch \"\$file\"; done"`)
    println("Submitted all batch jobs from: $dir_name")
end

function aggregate_first_stage_lb(first_lb_dir::String)
    println("Aggregating first-stage results...")
    
    # Find all simdirs
    simdirs = [joinpath(first_lb_dir, d) for d in readdir(first_lb_dir) 
               if isdir(joinpath(first_lb_dir, d)) && d != "batch_files"]
    
    if isempty(simdirs)
        error("No simulation directories found in $first_lb_dir")
    end
    
    # Read objectives from result files
    objectives = Float64[]
    
    for simdir in simdirs
        obj_file = joinpath(simdir, "output", "objective.csv")
        if isfile(obj_file)
            df = CSV.read(obj_file, DataFrame)
            push!(objectives, df[1, :objective])
        else
            error("Missing objective file: $obj_file. Job may not have completed yet.")
        end
    end
    
    # Get the largest objective
    first_lb_val = maximum(objectives)
    
    # Write it to a CSV
    df = DataFrame(first_stage_lb = first_lb_val)
    CSV.write(joinpath(first_lb_dir, "first_stage_lb.csv"), df)
    
    println("First-stage LB: $first_lb_val")
    return first_lb_val
    
end

function aggregate_second_stage_lb(second_lb_dir::String, superdir::String)
    println("Aggregating second-stage results...")
    
    # Read config to get representative probabilities
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    rep_prob = toml_data["representative_prob"]
    dates = toml_data["dates"]
    
    # Find simdirs in order matching dates
    simdirs = [joinpath(second_lb_dir, date) for date in dates]
    @assert length(simdirs) == length(rep_prob) "Mismatch between simdirs and probabilities"
    
    # Read objectives from result files
    objectives = Float64[]
    for simdir in simdirs
        obj_file = joinpath(simdir, "output", "objective.csv")
        if isfile(obj_file)
            df = CSV.read(obj_file, DataFrame)
            push!(objectives, df[1, :objective])
        else
            error("Missing objective file: $obj_file. Job may not have completed yet.")
        end
    end
    
    # Compute weighted average (order is preserved from dates)
    second_lb_val = sum(objectives[i] * rep_prob[i] for i in 1:length(simdirs))
    
    # Write to CSV
    df = DataFrame(second_stage_lb = second_lb_val)
    CSV.write(joinpath(second_lb_dir, "second_stage_lb.csv"), df)
    
    println("Second-stage LB: $second_lb_val")
    return second_lb_val
end

function compute_first_stage_lb(superdir, first_lb_dir; submit_jobs::Bool=true)
    """
    First-stage LB: computed by achieving zero load-shed feasibility on the "most challenging" representative day
    Solves <only_feasibility = true> for all representative days
    """
    # Read config file
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    
    # Create first stage simdirs
    simdirs = create_simdirs(first_lb_dir, toml_data, only_feasibility=true)
    
    # Create and submit batch files
    pace_dir = joinpath(first_lb_dir, "batch_files")
    create_batch_files(simdirs, pace_dir, "_first_lb", "first_stage", submit_jobs)
    
    println("First-stage batch jobs submitted. Results will be available after jobs complete.")
    println("To aggregate results later, run: aggregate_first_stage_lb(\"$first_lb_dir\")")
    
    return nothing
end

function compute_second_stage_lb(superdir, second_lb_dir; submit_jobs::Bool=true)
    """
    Second-stage LB: computed by taking the maximum parameters of all first-stage variables and optimizing operational cost
    Returns the weighted average of operational costs across representative days
    """
    # Read config file
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    
    gamma_max = toml_data["num_cap_upgrades_max"]
    s_energy_max = toml_data["max_energy_rating"]
    rep_prob = toml_data["representative_prob"]
    
    # Create second stage simdirs (order matches rep_prob)
    simdirs = create_simdirs(second_lb_dir, toml_data)
    @assert length(simdirs) == length(rep_prob) "Mismatch between simdirs and probabilities"
    
    # Write second stage parameters to config files
    for simdir in simdirs
        config_file = joinpath(simdir, "config.toml")
        config = TOML.parsefile(config_file)
        config["gamma_max"] = gamma_max
        config["s_energy_max"] = s_energy_max
        open(config_file, "w") do io
            TOML.print(io, config)
        end
    end
    
    # Create and submit batch files
    pace_dir = joinpath(second_lb_dir, "batch_files")
    create_batch_files(simdirs, pace_dir, "_second_lb", "second_stage", submit_jobs)
    
    println("Second-stage batch jobs submitted. Results will be available after jobs complete.")
    println("To aggregate results later, run: aggregate_second_stage_lb(\"$second_lb_dir\", \"$superdir\")")
    
    return nothing
end

function compute_benders_lb(superdir::String; force_lb::Bool=false)
    """
    Computes a valid global lower bound by summing first-stage and second-stage LBs.
    If CSV files exist, reads and returns their sum. Otherwise, submits batch jobs.
    
    # Arguments
    - `superdir::String`: Directory containing the main config file
    - `force_lb::Bool=false`: Force recomputation even if results exist
    
    # Returns
    - `Float64`: Total lower bound (first_lb + second_lb), or `nothing` if jobs need to run
    """
    lb_dir = mkpath(joinpath(superdir, "lower_bound"))
    first_lb_dir = mkpath(joinpath(lb_dir, "first_stage"))
    second_lb_dir = mkpath(joinpath(lb_dir, "second_stage"))
    first_lb_file = joinpath(first_lb_dir, "first_stage_lb.csv")
    second_lb_file = joinpath(second_lb_dir, "second_stage_lb.csv")
    
    # Check if both files exist and we're not forcing recomputation
    if !force_lb && isfile(first_lb_file) && isfile(second_lb_file)
        println("="^80)
        println("READING EXISTING LOWER BOUND RESULTS")
        println("="^80)
        
        # Read first-stage lower bound
        first_lb_df = CSV.read(first_lb_file, DataFrame)
        first_lb = first_lb_df[1, :first_stage_lb]  # Adjust column name as needed
        println("First-stage LB: $(round(first_lb, digits=2))")
        
        # Read second-stage lower bound
        second_lb_df = CSV.read(second_lb_file, DataFrame)
        second_lb = second_lb_df[1, :second_stage_lb]  # Adjust column name as needed
        println("Second-stage LB: $(round(second_lb, digits=2))")
        
        # Compute total
        total_lb = first_lb + second_lb
        println("\nTotal Benders LB: $(round(total_lb, digits=2))")
        println("="^80)
        
        return total_lb
    end
    
    # Otherwise, submit batch jobs
    println("="^80)
    println("BATCH MODE: Submitting lower bound computation jobs")
    println("="^80)
    
    # Submit first-stage jobs if needed
    if force_lb || !isfile(first_lb_file)
        println("\nSubmitting first-stage lower bound jobs...")
        compute_first_stage_lb(superdir, first_lb_dir)
    else
        println("\nFirst-stage results already exist at: $first_lb_file")
        println("Use force_lb=true to recompute")
    end
    
    # Submit second-stage jobs if needed
    if force_lb || !isfile(second_lb_file)
        println("\nSubmitting second-stage lower bound jobs...")
        compute_second_stage_lb(superdir, second_lb_dir)
    else
        println("\nSecond-stage results already exist at: $second_lb_file")
        println("Use force_lb=true to recompute")
    end
    
    println("\n" * "="^80)
    println("BATCH JOBS SUBMITTED")
    println("="^80)
    println("\nAfter jobs complete, run this function again to compute the bound.")
    println("Or manually aggregate with:")
    println("  first_lb = aggregate_first_stage_lb(\"$first_lb_dir\")")
    println("  second_lb = aggregate_second_stage_lb(\"$second_lb_dir\", \"$superdir\")")
    println("  benders_lb = first_lb + second_lb")
    println("="^80)
    
    return nothing
end