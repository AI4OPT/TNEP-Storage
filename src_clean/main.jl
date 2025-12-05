using JuMP
include("structs/ExpansionPlanner.jl")
include("structs/ParallelizedBenders.jl")
include("utils/utils.jl")

function parallelized_benders_wrapper(superdir::String; 
                                  max_iterations::Int=100000,
                                  tolerance::Float64=0.01,
                                  force_lb::Bool=false)
    println("[DEBUG] Computing global lower bound...")
    lb_value = compute_benders_lb(superdir, force_lb=force_lb)

    # benders = parallelized_ptdf_benders(superdir, max_iterations=max_iterations, tolerance=tolerance, )
end

# superdir = "sim/PowerUp/benders/2030batchtest"
# lb_value = compute_benders_lb(superdir)
# parallelized_ptdf_benders(superdir, lower_bound=lb_value)