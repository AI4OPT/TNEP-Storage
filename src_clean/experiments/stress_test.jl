"""
    setup_stress_test(stressdir; 
                      pace_dir=nothing,
                      name_suffix="",
                      submit_jobs=false)

Create simulation directories and optionally generate SLURM batch files for stress testing.

# Arguments
- `stressdir`: Base directory containing config.toml
- `pace_dir`: Directory for SLURM batch files (optional, required if generating batch files)
- `name_suffix`: Suffix to append to job names (default: "")
- `submit_jobs`: Whether to automatically submit generated batch jobs (default: false)

# Returns
- Vector of created simulation directory paths
"""
function setup_stress_test(stressdir::String;
                          pace_dir::Union{String,Nothing}=nothing,
                          name_suffix::String="",
                          submit_jobs::Bool=false)
    # Load and validate config
    config_file = joinpath(stressdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    
    # Check for investment files
    investment_files = [
        "line_investments.csv",
        "storage_investments.csv"
    ]
    available_investments = filter(f -> isfile(joinpath(stressdir, f)), investment_files)
    
    # Create simulation directories
    simdirs = create_simulation_directories(stressdir, toml_data, available_investments)
    
    # Optionally create batch files
    if !isnothing(pace_dir)
        create_batch_files(simdirs, pace_dir, name_suffix, submit_jobs)
    end
    
    return simdirs
end

"""
Create individual simulation directories for each date in the config.
"""
function create_simulation_directories(stressdir::String, 
                                      toml_data::Dict,
                                      investment_files::Vector{String})
    simdirs = String[]
    
    for rep in toml_data["dates"]
        simdir = joinpath(stressdir, rep)
        
        # Create directory if needed
        mkpath(simdir)
        
        # Create modified config
        config = deepcopy(toml_data)
        config["dates"] = [rep]
        config["representative_prob"] = [1.0]
        config["num_representatives"] = 1
        
        # Copy investment files
        for file in investment_files
            src = joinpath(stressdir, file)
            dst = joinpath(simdir, file)
            cp(src, dst; force=true)
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

"""
Create SLURM batch files for all simulation directories.
"""
function create_batch_files(simdirs::Vector{String},
                           pace_dir::String,
                           name_suffix::String,
                           submit_jobs::Bool)
    mkpath(pace_dir)
    
    for simdir in simdirs
        rep = basename(simdir)
        job_name = "$(rep)$(name_suffix)"
        log_file = joinpath("stress_test_2", job_name)
        output_file = joinpath(pace_dir, "run_$(rep).sbatch")
        
        create_sbatch_file(job_name, log_file, simdir, output_file)
    end
    
    if submit_jobs
        submit_sbatch_jobs(pace_dir)
    end
end

"""
Generate a SLURM batch file with the specified parameters.
"""
function create_sbatch_file(job_name::String, 
                           log_file::String, 
                           sim_path::String, 
                           output_file::String;
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
    julia --project=. exp/clean/run_expansionplanner.jl $sim_path $tasks_per_node
    """
    
    open(output_file, "w") do file
        write(file, content)
    end
    
    println("Created sbatch file: $output_file")
end

"""
Submit all .sbatch files in the specified directory.
"""
function submit_sbatch_jobs(dir_name::String)
    run(`bash -c "for file in \"$dir_name\"/*.sbatch; do sbatch \"\$file\"; done"`)
    println("Submitted all batch jobs from: $dir_name")
end

# Usage examples:
# 1. Just create simulation directories
# setup_stress_test("path/to/stress/dir")

# 2. Create directories and batch files (don't submit)
# setup_stress_test("path/to/stress/dir"; 
#                   pace_dir="PACE/r1/PowerUp/2045",
#                   name_suffix="_2045")

# 3. Create directories, batch files, and submit jobs
# setup_stress_test("path/to/stress/dir"; 
#                   pace_dir="PACE/r1/PowerUp/2045",
#                   name_suffix="_2045",
#                   submit_jobs=true)

"""
stressdir
pace_dir = "PACE/r1/PowerUp/2030nobenders"
name_suffix="_2030nobenders"

setup_stress_test(stressdir, pace_dir=pace_dir, name_suffix=name_suffix, submit_jobs=true)
"""