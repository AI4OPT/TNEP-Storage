using TOML
using CSV
using DataFrames

include("../structs/ExpansionPlanner.jl")
most_challenging = "2016-08-11"

function create_sbatch_file(job_name::String, 
                           log_file::String, 
                           sim_path::String, 
                           output_file::String,
                           nodes::Int=1,
                           tasks_per_node::Int=24,
                           mem_per_cpu::String="12G",
                           time_limit::String="24:00:00",
                           partition::String="inferno",
                           account::String="gts-phentenryck3-ai4opt",
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
    julia --project=. exp/clean/run_lb_single.jl $sim_path $tasks_per_node
    """
    open(output_file, "w") do file
        write(file, content)
    end
    println("Created sbatch file: $output_file")
end

function create_batch_files(simdirs::Vector{String},
                           pace_dir::String,
                           name_suffix::String,
                           submit_jobs::Bool)
    mkpath(pace_dir)
    for simdir in simdirs
        rep = basename(simdir)
        job_name = "$(rep)$(name_suffix)"
        log_file = joinpath("lower_bound", job_name)
        output_file = joinpath(pace_dir, "run_$(rep).sbatch")
        create_sbatch_file(job_name, log_file, simdir, output_file)
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
        summary = joinpath(simdir, "output", "summary_data.csv")
        if isfile(summary)
            df = CSV.read(summary, DataFrame)
        
            # Extract line and storage investment costs
            line_cost_row = df[df.Variable .== "line_investment_costs", :]
            storage_cost_row = df[df.Variable .== "storage_investment_costs", :]
        
            line_cost = isempty(line_cost_row) ? 0.0 : line_cost_row[1, :Value]
            storage_cost = isempty(storage_cost_row) ? 0.0 : storage_cost_row[1, :Value]
        
            # Push the sum of both investment costs
            push!(objectives, line_cost + storage_cost)
        else
            error("Missing objective file: $summary. Job may not have completed yet.")
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
    year = toml_data["decarbonization_year"]
    real_dates = [string(year) * date[5:end] for date in dates]
    
    # Find simdirs in order matching dates
    simdirs = [joinpath(second_lb_dir, date) for date in real_dates]
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

function compute_first_stage_lb(superdir, first_lb_dir, most_challenging; submit_jobs::Bool=true)
    """
    First-stage LB: computed by achieving zero load-shed feasibility on the "most challenging" representative day
    E.g. "most challenging" == "2016-08-11
    Solves <only_feasibility = true> the most challenging representative day
    """
    # Read config file
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    
    # Create the simdir
    year = toml_data["decarbonization_year"]
    simdir = joinpath(first_lb_dir, string(year) * most_challenging[5:end])
    mkpath(simdir)
    
    config = deepcopy(toml_data)
    config["dates"] = [most_challenging]
    config["representative_prob"] = [1.0]
    config["num_representatives"] = 1
    config["only_feasibility"] = true
    
    # Write config
    open(joinpath(simdir, "config.toml"), "w") do io
        TOML.print(io, config)
    end
    
    # Create batch file
    pace_dir = joinpath(first_lb_dir, "batch_files")
    mkpath(pace_dir)
    
    job_name = "$(basename(simdir))_first_lb"
    log_file = joinpath("lower_bound", job_name)
    batch_file = joinpath(pace_dir, "$(job_name).sbatch")
    
    content = """
    #!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-ai4opt
    #SBATCH -N1 --ntasks-per-node=24
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t48:00:00
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$log_file.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu
    
    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=. exp/clean/run_expansionplanner.jl $simdir
    """
    
    open(batch_file, "w") do file
        write(file, content)
    end
    
    # Submit the job if requested
    if submit_jobs
        run(`sbatch $batch_file`)
        println("First-stage batch job submitted: $batch_file")
    else
        println("Batch file created but not submitted: $batch_file")
    end
    
    println("Results will be available after job completes.")
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
    create_batch_files(simdirs, pace_dir, "_second_lb", submit_jobs)
    
    println("Second-stage batch jobs submitted. Results will be available after jobs complete.")
    println("To aggregate results later, run: aggregate_second_stage_lb(\"$second_lb_dir\", \"$superdir\")")
    
    return nothing
end

function compute_benders_lb(superdir::String, most_challenging::String; force_lb::Bool=false, submit_jobs::Bool=true)
    """
    Computes a valid global lower bound by summing first-stage and second-stage LBs.
    Priority order:
    1. If CSV files exist, reads and returns their sum
    2. If result files exist but not aggregated, aggregates them
    3. Otherwise, submits batch jobs (if submit_jobs=true)
    
    # Arguments
    - `superdir::String`: Directory containing the main config file
    - `most_challenging::String`: Most challenging representative day
    - `force_lb::Bool=false`: Force recomputation even if results exist
    - `submit_jobs::Bool=true`: Whether to submit batch jobs if needed
    
    # Returns
    - `Float64`: Total lower bound (first_lb + second_lb), or `nothing` if jobs need to run
    """
    lb_dir = mkpath(joinpath(superdir, "lower_bound"))
    first_lb_dir = mkpath(joinpath(lb_dir, "first_stage"))
    second_lb_dir = mkpath(joinpath(lb_dir, "second_stage"))
    first_lb_file = joinpath(first_lb_dir, "first_stage_lb.csv")
    second_lb_file = joinpath(second_lb_dir, "second_stage_lb.csv")
    
    # Priority 1: Check if both aggregated files exist and we're not forcing recomputation
    if !force_lb && isfile(first_lb_file) && isfile(second_lb_file)
        println("="^80)
        println("READING EXISTING LOWER BOUND RESULTS")
        println("="^80)
        
        # Read first-stage lower bound
        first_lb_df = CSV.read(first_lb_file, DataFrame)
        first_lb = first_lb_df[1, :first_stage_lb]
        println("First-stage LB: $(round(first_lb, digits=2))")
        
        # Read second-stage lower bound
        second_lb_df = CSV.read(second_lb_file, DataFrame)
        second_lb = second_lb_df[1, :second_stage_lb]
        println("Second-stage LB: $(round(second_lb, digits=2))")
        
        # Compute total
        total_lb = first_lb + second_lb
        println("\nTotal Benders LB: $(round(total_lb, digits=2))")
        println("="^80)
        
        return total_lb
    end
    
    # Priority 2: Try to aggregate if result files exist
    first_lb = nothing
    second_lb = nothing
    
    if !force_lb
        # Check if first-stage results can be aggregated
        if !isfile(first_lb_file)
            year = TOML.parsefile(joinpath(superdir, "config.toml"))["decarbonization_year"]
            first_stage_simdir = joinpath(first_lb_dir, string(year) * most_challenging[5:end])
            summary_file = joinpath(first_stage_simdir, "output", "summary_data.csv")
            
            if isfile(summary_file)
                println("="^80)
                println("AGGREGATING FIRST-STAGE RESULTS")
                println("="^80)
                try
                    first_lb = aggregate_first_stage_lb(first_lb_dir)
                catch e
                    println("Warning: Could not aggregate first-stage results: $e")
                end
            end
        else
            first_lb_df = CSV.read(first_lb_file, DataFrame)
            first_lb = first_lb_df[1, :first_stage_lb]
            println("First-stage LB already computed: $(round(first_lb, digits=2))")
        end
        
        # Check if second-stage results can be aggregated
        if !isfile(second_lb_file)
            config = TOML.parsefile(joinpath(superdir, "config.toml"))
            year = config["decarbonization_year"]
            dates = config["dates"]
            real_dates = [string(year) * date[5:end] for date in dates]
            
            # Check if at least one objective file exists
            any_obj_exists = any(isfile(joinpath(second_lb_dir, date, "output", "objective.csv")) for date in real_dates)
            
            if any_obj_exists
                println("="^80)
                println("AGGREGATING SECOND-STAGE RESULTS")
                println("="^80)
                try
                    second_lb = aggregate_second_stage_lb(second_lb_dir, superdir)
                catch e
                    println("Warning: Could not aggregate second-stage results: $e")
                end
            end
        else
            second_lb_df = CSV.read(second_lb_file, DataFrame)
            second_lb = second_lb_df[1, :second_stage_lb]
            println("Second-stage LB already computed: $(round(second_lb, digits=2))")
        end
        
        # If both bounds were successfully obtained, return the sum
        if first_lb !== nothing && second_lb !== nothing
            total_lb = first_lb + second_lb
            println("\n" * "="^80)
            println("Total Benders LB: $(round(total_lb, digits=2))")
            println("="^80)
            return total_lb
        end
    end
    
    # Priority 3: Submit batch jobs if enabled
    if !submit_jobs
        println("="^80)
        println("CANNOT COMPUTE LOWER BOUND")
        println("="^80)
        println("Results not available and submit_jobs=false")
        println("Set submit_jobs=true to submit batch jobs")
        return nothing
    end
    
    println("="^80)
    println("BATCH MODE: Submitting lower bound computation jobs")
    println("="^80)
    
    # Submit first-stage jobs if needed
    if force_lb || first_lb === nothing
        println("\nSubmitting first-stage lower bound jobs...")
        compute_first_stage_lb(superdir, first_lb_dir, most_challenging, submit_jobs=true)
    end
    
    # Submit second-stage jobs if needed
    if force_lb || second_lb === nothing
        println("\nSubmitting second-stage lower bound jobs...")
        compute_second_stage_lb(superdir, second_lb_dir, submit_jobs=true)
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