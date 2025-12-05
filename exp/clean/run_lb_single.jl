using TOML
using CSV
using DataFrames

include("../../src_clean/main.jl")

"""
Run a single simulation for lower bound computation.
This script is called by SLURM batch jobs.

Usage: julia run_lb_single.jl <sim_path> <num_threads> <stage_type>
"""

if length(ARGS) < 3
    error("Usage: julia run_lb_single.jl <sim_path> <num_threads> <stage_type>")
end

sim_path = ARGS[1]
num_threads = parse(Int, ARGS[2])
stage_type = ARGS[3]  # "first_stage" or "second_stage"

println("="^80)
println("Running lower bound computation:")
println("  Simulation path: $sim_path")
println("  Threads: $num_threads")
println("  Stage type: $stage_type")
println("="^80)

# Read config
config_file = joinpath(sim_path, "config.toml")
!isfile(config_file) && error("Config file not found: $config_file")
config = TOML.parsefile(config_file)

if stage_type == "first_stage"
    # First stage: run with only_feasibility=true
    println("Running first-stage computation (feasibility)...")
    
    # The config should already have only_feasibility set, but verify
    if !haskey(config, "only_feasibility") || !config["only_feasibility"]
        println("Warning: only_feasibility not set in config, setting it now")
        config["only_feasibility"] = true
        open(config_file, "w") do io
            TOML.print(io, config)
        end
    end
    
    # Create planner and run model
    planner = ExpansionPlanner(sim_path)
    create_model!(planner)
    configure_optimizer!(planner)
    solve!(planner.model)
    
    # Check status
    jump_model = planner.model.jump_model
    if primal_status(jump_model) != MOI.FEASIBLE_POINT
        error("Infeasible solution for $sim_path")
    end
    if termination_status(jump_model) != MOI.OPTIMAL
        println("Warning: Termination status not optimal")
    end
    
    # Export results
    export_results!(planner)
    
    # Export objective value
    obj_val = objective_value(jump_model)
    results_dir = mkpath(joinpath(sim_path, "output"))
    df = DataFrame(objective = obj_val)
    CSV.write(joinpath(results_dir, "objective.csv"), df)
    println("First-stage objective: $obj_val")
    
elseif stage_type == "second_stage"
    # Second stage: fix investments to max and optimize operations
    println("Running second-stage computation (operational cost)...")
    
    # Get parameters from config
    gamma_max = config["gamma_max"]
    s_energy_max = config["s_energy_max"]
    
    # Create and solve model
    planner = ExpansionPlanner(sim_path)
    create_model!(planner)
    configure_optimizer!(planner)
    
    # Fix maximum parameters of all first-stage variables
    gamma_val = fill(gamma_max, planner.model.E)
    s_energy_val = fill(s_energy_max, planner.model.N)
    fix_investments!(planner.model, gamma_val, s_energy_val)
    set_objective!(planner.model, :operational)
    solve!(planner.model)
    
    # Check status
    jump_model = planner.model.jump_model
    if primal_status(jump_model) != MOI.FEASIBLE_POINT
        error("Infeasible solution for $sim_path")
    end
    if termination_status(jump_model) != MOI.OPTIMAL
        println("Warning: Termination status not optimal")
    end
    
    # Export results
    export_results!(planner)
    
    # Export objective value
    obj_val = objective_value(jump_model)
    results_dir = mkpath(joinpath(sim_path, "output"))
    df = DataFrame(objective = obj_val)
    CSV.write(joinpath(results_dir, "objective.csv"), df)
    println("Second-stage objective: $obj_val")
    
else
    error("Unknown stage_type: $stage_type. Must be 'first_stage' or 'second_stage'")
end

println("="^80)
println("Computation completed successfully!")
println("="^80)