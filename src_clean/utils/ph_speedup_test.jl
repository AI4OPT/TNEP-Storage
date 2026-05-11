using Gurobi, JuMP, CSV, DataFrames

include("../structs/ExpansionPlanner.jl")
include("../structs/PHSubproblem.jl")

const PROJECT_DIR = "/storage/home/hcoda1/1/kwu381/TNEP-Storage"

# ---------------------------------------------------------------------------
# run_ph_speedup_test
#
# Creates and submits a single SBATCH job for a given simdir and mode:
#   mode="full" — monolithic PTDFModel solve via run_model (uses
#                 current_investment_dir from params to fix investments)
#   mode="ph"   — PH temporal decomposition
#
# Results are written to simdir/ph_speedup/result_{mode}.csv.
# ---------------------------------------------------------------------------
function run_ph_speedup_test(simdir::String, mode::String;
                              submit_jobs::Bool=true,
                              n_blocks::Int=-1,
                              max_ph_iter::Int=-1,
                              ph_tol::Float64=-1.0)

    mode ∈ ("full", "ph") || error("mode must be 'full' or 'ph', got '$mode'")

    _create_speedup_sbatch(simdir, mode;
                           submit_jobs=submit_jobs,
                           n_blocks=n_blocks,
                           max_ph_iter=max_ph_iter, ph_tol=ph_tol)
end

function _create_speedup_sbatch(simdir::String, mode::String;
                                 submit_jobs::Bool=true,
                                 n_blocks::Int=-1,
                                 max_ph_iter::Int=-1,
                                 ph_tol::Float64=-1.0)
    job_name = basename(simdir) * "_ph_speedup_" * mode
    log_dir  = mkpath(joinpath(PROJECT_DIR, "PACE", "logs", "ph_speedup"))
    script   = joinpath(PROJECT_DIR, "src_clean", "utils", "ph_speedup_test.jl")

    content = """
    #!/bin/bash
    #SBATCH -J$job_name
    #SBATCH -qinferno
    #SBATCH --account=gts-phentenryck3-ai4opt
    #SBATCH -N1 --ntasks-per-node=4
    #SBATCH --mem-per-cpu=16G
    #SBATCH -t12:00:00
    #SBATCH -o$(joinpath(log_dir, job_name)).out
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=kwu381@gatech.edu

    cd $PROJECT_DIR
    julia --project=. $script $simdir $mode $n_blocks $max_ph_iter $ph_tol
    """

    batch_dir   = mkpath(joinpath(simdir, "batch_file"))
    output_file = joinpath(batch_dir, "$job_name.sbatch")
    open(output_file, "w") do f; write(f, content); end
    println("Created: $output_file")

    if submit_jobs
        run(`sbatch $output_file`)
        println("Submitted: $job_name")
    end

    return output_file
end

# ---------------------------------------------------------------------------
# Script entry point — called by SBATCH with positional args:
#   julia ph_speedup_test.jl <simdir> <mode> [n_blocks] [rho] [max_ph_iter] [ph_tol]
#
# For mode="full": current_investment_dir must be set in the simdir params.json
#   so that PTDFModel auto-fixes investments (see PTDFModel constructor).
# ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) < 2 && error("Usage: julia ph_speedup_test.jl <simdir> <mode> [n_blocks max_ph_iter ph_tol]")

    simdir      = ARGS[1]
    mode        = ARGS[2]
    n_blocks    = length(ARGS) >= 3 ? parse(Int,     ARGS[3]) : -1
    max_ph_iter = length(ARGS) >= 4 ? parse(Int,     ARGS[4]) : -1
    ph_tol      = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : -1.0

    output_dir = mkpath(joinpath(simdir, "ph_speedup"))

    if mode == "full"
        println("\n=== Full PTDFModel (via run_model) ===")
        t0 = time()
        jump_model, _ = run_model(simdir)
        t = time() - t0
        obj = objective_value(jump_model)

        tracked_csv = joinpath(simdir, "output", "tracked_constraints.csv")
        n_cuts = isfile(tracked_csv) ? nrow(CSV.read(tracked_csv, DataFrame)) : 0

        println("time=$(round(t, digits=2))s | obj=$(round(obj, digits=2)) | ptdf_cuts=$n_cuts")
        CSV.write(joinpath(output_dir, "result_full.csv"), DataFrame(
            mode="full", solve_time=t, obj=obj,
            n_ptdf_cuts=n_cuts, n_blocks=1, rho=0.0, boundary_residual=0.0
        ))

    elseif mode == "ph"
        data    = set_up_data(simdir)
        inv_dir = get(data["param"], "current_investment_dir", joinpath(simdir, "output"))
        E, N    = length(data["branch"]), length(data["bus"])
        gamma_val, s_energy_val = load_investments_from_dir(inv_dir, data; E=E, N=N, allow_missing=false)

        println("Loaded investments: $(sum(gamma_val .> 0)) upgraded lines, $(sum(s_energy_val .> 0)) storage nodes")

        ph = PHSubproblem(data, Gurobi.Optimizer, simdir;
                                  n_blocks=n_blocks,
                                  max_ph_iterations=max_ph_iter, ph_tol=ph_tol)
        println("\n=== PH temporal decomposition (n_blocks=$(ph.n_blocks), ph_rho=$(get(data["param"],"ph_rho",1.0)), ph_obj_scale=$(get(data["param"],"ph_obj_scale",1.0))) ===")
        fix_investments!(ph, gamma_val, s_energy_val)

        t0  = time()
        solve!(ph)
        t   = time() - t0
        obj = get_objective_value(ph)
        n_cuts   = sum(length(ph.block_models[r,b].tracked_constraints)
                       for r in 1:ph.R, b in 1:ph.n_blocks)
        residual = maximum(abs.(ph.soc_out .- ph.soc_in))

        println("time=$(round(t, digits=2))s | obj=$(round(obj, digits=2)) | ptdf_cuts=$n_cuts | residual=$(round(residual, digits=6))")
        CSV.write(joinpath(output_dir, "result_ph.csv"), DataFrame(
            mode="ph", solve_time=t, obj=obj,
            n_ptdf_cuts=n_cuts, n_blocks=ph.n_blocks,
            rho=get(data["param"], "ph_rho", 1.0),
            obj_scale=get(data["param"], "ph_obj_scale", 1.0),
            boundary_residual=residual
        ))

    else
        error("Unknown mode '$mode' — use 'full' or 'ph'")
    end
end
