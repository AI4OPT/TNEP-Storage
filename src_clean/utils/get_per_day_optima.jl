using TOML
using CSV
using DataFrames

include("../structs/ExpansionPlanner.jl")

function get_per_day_optima(superdir::String; submit_jobs::Bool=true)
    # Load and validate config
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)

    simdirs = create_simdirs(superdir, toml_data)

    for simdir in simdirs
        create_per_day_sbatch_file(simdir, submit_jobs=submit_jobs)
    end
end

function create_per_day_sbatch_file(simdir::String; submit_jobs::Bool=true)
    job_name = basename(simdir) * "_per_day"
    log_file = joinpath("per_day", job_name)

    content = """
    #!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-coda20
    #SBATCH -N1 --ntasks-per-node=24
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t24:00:00
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$log_file.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu

    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=. exp/clean/run_expansionplanner.jl $simdir
    """

    batch_dir = mkpath(joinpath(simdir, "batch_file"))
    output_file = joinpath(batch_dir, "$job_name.sbatch")

    open(output_file, "w") do file
        write(file, content)
    end
    println("Created sbatch file: $output_file")

    if submit_jobs
        run(`sbatch $output_file`)
        println("Submitted job: $job_name")
    end

    return output_file
end