function write_parallelizedbenders_sbatch_file(simdir; threads=1, hours=24, submit=true)
    pace_dir = "PACE/r1/PowerUp/parallelizedbenders"
    job_name = basename(simdir)

    n_nodes = 1
    if threads > 24
        threads = 24
    end
    
    # Create the SBATCH content
    sbatch_content = """#!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-coda20
    #SBATCH -N$n_nodes --ntasks-per-node=$threads
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t$hours:00:00
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$job_name.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu

    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=. -t $threads exp/clean/run_parallelizedbenders.jl $simdir
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

function write_expansionplanner_sbatch_file(simdir; hours=24, submit=true)
    pace_dir = "PACE/r1/PowerUp/expansionplanner"
    job_name = basename(simdir)
    
    # Create the SBATCH content
    sbatch_content = """#!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-coda20
    #SBATCH -N1 --ntasks-per-node=24
    #SBATCH --mem-per-cpu=12G
    #SBATCH -t$hours:00:00
    #SBATCH -o/storage/home/hcoda1/1/kwu381/TNEP-Storage/PACE/logs/$job_name.out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu

    cd /storage/home/hcoda1/1/kwu381/TNEP-Storage
    julia --project=.  exp/clean/run_expansionplanner.jl $simdir
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