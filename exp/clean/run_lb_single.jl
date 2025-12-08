using TOML
using CSV
using DataFrames
include("../../src_clean/main.jl")

"""
Run second-stage lower bound computation with fixed maximum investments.
This script is called by SLURM batch jobs.
Usage: julia run_lb_single.jl <sim_path> <num_threads>
"""

if length(ARGS) < 2
    error("Usage: julia run_lb_single.jl <sim_path> <num_threads>")
end

sim_path = ARGS[1]
num_threads = parse(Int, ARGS[2])

println("="^80)
println("Running second-stage lower bound computation (operational cost):")
println("  Simulation path: $sim_path")
println("  Threads: $num_threads")
println("="^80)

# Read config
config_file = joinpath(sim_path, "config.toml")
!isfile(config_file) && error("Config file not found: $config_file")
config = TOML.parsefile(config_file)

# Get parameters from config
gamma_max = config["gamma_max"]
s_energy_max = config["s_energy_max"]

# Create and solve model
planner = ExpansionPlanner(sim_path)
create_model!(planner)
configure_optimizer!(planner)

# Fix investments to maximum values
gamma_val = fill(gamma_max, planner.model.E)
s_energy_val = fill(s_energy_max, planner.model.N)
fix_investments!(planner.model, gamma_val, s_energy_val)

# Set operational objective
set_objective!(planner.model, :operational)
solve!(planner.model)

# Check solution status
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
println("="^80)
println("Computation completed successfully!")
println("="^80)