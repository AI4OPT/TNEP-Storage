# function write_parallelizedbenders_sbatch_file(simdir; threads=1, hours=24, submit=true)
#     pace_dir = "PACE/r1/PowerUp2026Rebuttal/parallelizedbenders"
#     job_name = basename(simdir)

#     n_nodes = 1
#     if threads > 24
#         threads = 24
#     end
    
#     # Create the SBATCH content
#     sbatch_content = """#!/bin/bash
#     #SBATCH -J$job_name
#     #SBATCH -qinferno
#     #SBATCH --account=gts-phentenryck3-ai4opt
#     #SBATCH -N$n_nodes --ntasks-per-node=$threads
#     #SBATCH --mem-per-cpu=12G
#     #SBATCH -t$hours:00:00
#     #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$job_name.out
#     #SBATCH --mail-type=BEGIN,END,FAIL
#     #SBATCH --mail-user=kwu381@gatech.edu

#     cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
#     julia --project=. -t $threads exp/clean/run_parallelizedbenders.jl $simdir
#     """
    
#     # Write to file
#     sbatch_file = joinpath(pace_dir, "$job_name.sbatch")

#     open(sbatch_file, "w") do file
#         write(file, sbatch_content)
#     end
    
#     println("SBATCH file written to: $sbatch_file")

#     if submit
#         run(`sbatch $sbatch_file`)
#     end

#     return sbatch_file
# end

function write_distributed_parallelizedbenders_sbatch_file(simdir; n_workers=nothing, hours=24, submit=true)
    pace_dir = "PACE/r1/PowerUp2026Rebuttal/distributedparallelizedbenders"
    job_name = basename(simdir)
    
    # Determine number of workers (should equal number of rep days)
    if n_workers === nothing
        # Auto-detect from config
        config_file = joinpath(simdir, "config.toml")
        toml_data = TOML.parsefile(config_file)
        n_workers = length(toml_data["dates"]) * length(get(toml_data, "years", ["0000"])) + 1

    end
    
    # SLURM configuration
    max_tasks_per_node = 24
    
    # Calculate nodes needed (ceiling division)
    n_nodes = cld(n_workers, max_tasks_per_node)
    
    # Calculate tasks per node (distribute evenly)
    tasks_per_node = cld(n_workers, n_nodes)
    
    # Create the SBATCH content
    sbatch_content = """#!/bin/bash
    #SBATCH -J $job_name
    #SBATCH -q inferno
    #SBATCH --account=gts-phentenryck3-ai4opt
    #SBATCH -N $n_nodes
    #SBATCH --ntasks=$n_workers
    #SBATCH --ntasks-per-node=$tasks_per_node
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t $hours:00:00
    #SBATCH -o /storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$job_name.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu

    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage

    # Run with use_slurm=true (will auto-detect SLURM environment)
    julia --project=. exp/clean/run_parallelizedbendersdistributed.jl $simdir
    """

    # Write to file
    mkpath(pace_dir)  # Ensure directory exists
    sbatch_file = joinpath(pace_dir, "$job_name.sbatch")
    
    open(sbatch_file, "w") do file
        write(file, sbatch_content)
    end
    
    println("SBATCH file written to: $sbatch_file")
    println("  Nodes: $n_nodes")
    println("  Total tasks: $n_workers")
    println("  Tasks per node: $tasks_per_node")
    
    if submit
        run(`sbatch $sbatch_file`)
    end
    
    return sbatch_file
end

function write_expansionplanner_sbatch_file(simdir; hours=24, submit=true, duals=false)
    pace_dir = "PACE/r1/PowerUp2026Rebuttal/expansionplanner"
    job_name = basename(simdir)
    
    duals_flag = duals ? " -d" : ""
    
    # Create the SBATCH content
    sbatch_content = """#!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-ai4opt
    #SBATCH -N1 --ntasks-per-node=24
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t$hours:00:00
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$job_name.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu
    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=.  exp/clean/run_expansionplanner.jl $simdir$(duals_flag)
    """
    # Write to file
    sbatch_file = joinpath(pace_dir, "$job_name.sbatch")
    open(sbatch_file, "w") do file
        write(file, sbatch_content)
    end
    println("SBATCH file written to: $sbatch_file")
    if submit
        run(`sbatch $sbatch_file`)
    end
    return sbatch_file
end