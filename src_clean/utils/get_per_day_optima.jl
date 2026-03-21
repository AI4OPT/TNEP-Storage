using TOML
using CSV
using DataFrames

include("../structs/ExpansionPlanner.jl")
include("../helpers/process_max_upgrades.jl")
include("../helpers/export_model.jl")

function get_per_day_optima(superdir::String; submit_jobs::Bool=true, force::Bool=false)
    # Load and validate config
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)

    simdirs = create_simdirs(superdir, toml_data)
    n_submitted_jobs = 0

    for simdir in simdirs
        if !isfile(joinpath(simdir, "output", "summary_data.csv")) || force
            create_per_day_sbatch_file(simdir, submit_jobs=submit_jobs)
            n_submitted_jobs += 1
        end
    end

    if n_submitted_jobs == 0
        compute_avg_summary(superdir, toml_data)
        process_max_aggregate(superdir, simdirs[1])
    end
    
end

function compute_avg_summary(superdir, toml_data)
    mkpath(joinpath(superdir, "output"))
    year = toml_data["decarbonization_year"]
    probs = toml_data["representative_prob"]
    weighted_sums = Dict{String, Float64}()
    
    for (i, rep) in enumerate(toml_data["dates"])
        simdir = joinpath(superdir, string(year) * rep[5:end])
        filepath = joinpath(simdir, "output", "summary_data.csv")
        df = CSV.read(filepath, DataFrame)
        
        # Get the weight for this representative day
        weight = probs[i]
        
        # Accumulate weighted values
        for row in eachrow(df)
            var_name = row.Variable
            value = row.Value
            
            # Skip NaN values
            if !isnan(value)
                if haskey(weighted_sums, var_name)
                    weighted_sums[var_name] += value * weight
                else
                    weighted_sums[var_name] = value * weight
                end
            end
        end
    end
    
    avg_df = DataFrame(
        Variable = collect(keys(weighted_sums)),
        Value = collect(values(weighted_sums))
    )
    sort!(avg_df, :Variable)
    
    CSV.write(joinpath(superdir, "output", "summary.csv"), avg_df)
    
    # Print sum of all cost variables
    total_costs = sum(row.Value for row in eachrow(avg_df) if endswith(row.Variable, "costs"))
    println("Total costs: $total_costs")
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

function process_max_aggregate(superdir, simdir)
    gamma_val, s_energy_val = compute_superset_core_point(superdir, is_multistage=false)
    data = JSON.parsefile(joinpath(simdir, "data.json"))
    export_investments_csv(data, gamma_val, s_energy_val, output_dir=joinpath(superdir, "output"), file_suffix="")
end