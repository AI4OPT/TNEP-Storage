using Base.Threads
using TOML
using CSV
using DataFrames

include("../structs/ExpansionPlanner.jl")

"""
Create individual simulation directories for each date in the config.
"""
function create_simdirs(superdir::String, toml_data::Dict; only_feasibility=nothing)
    simdirs = String[]
    
    for rep in toml_data["dates"]
        simdir = joinpath(superdir, rep)
        
        # Create directory if needed
        mkpath(simdir)
        
        # Create modified config
        config = deepcopy(toml_data)
        config["dates"] = [rep]
        config["representative_prob"] = [1.0]
        config["num_representatives"] = 1

        if !isnothing(only_feasibility)
            config["only_feasibility"] = only_feasibility
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

function compute_first_stage_lb(superdir, first_lb_dir)
    """
    First-stage LB: computed by achieving zero load-shed feasibility on the "most challenging" representative day
    Solves <only_feasibility = true> for all representative days
    """
    # Read config file
    config_file = joinpath(superdir, "config.toml")
    !isfile(config_file) && error("Config file not found: $config_file")
    toml_data = TOML.parsefile(config_file)
    
    # Create first stage simdirs
    simdirs = create_simdirs(first_lb_dir, toml_data, only_feasibility=true)
    
    # Solve all simdirs in parallel
    n_sims = length(simdirs)
    jump_models = Vector{Any}(undef, n_sims)
    datas = Vector{Any}(undef, n_sims)
    
    Threads.@threads for i in 1:n_sims
        simdir = simdirs[i]
        jump_models[i], datas[i] = run_model(simdir)
    end
    
    # Extract objectives from all solved models
    objectives = [objective_value(model) for model in jump_models]
    
    # Get the largest objective from all the jump_models
    first_lb_val = maximum(objectives)
    
    # Write it to a CSV
    df = DataFrame(first_stage_lb = first_lb_val)
    CSV.write(joinpath(first_lb_dir, "first_stage_lb.csv"), df)
    
    return first_lb_val, jump_models, datas
end

function compute_second_stage_lb(superdir, second_lb_dir)
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
    rep_prob = toml_data["representative_prob"]  # Weights in correct order
    
    # Create second stage simdirs (order matches rep_prob)
    simdirs = create_simdirs(second_lb_dir, toml_data)
    @assert length(simdirs) == length(rep_prob) "Mismatch between simdirs and probabilities"
    
    # Solve all simdirs in parallel
    n_sims = length(simdirs)
    jump_models = Vector{Any}(undef, n_sims)
    datas = Vector{Any}(undef, n_sims)
    objectives = Vector{Float64}(undef, n_sims)
    
    Threads.@threads for i in 1:n_sims
        simdir = simdirs[i]
        planner = ExpansionPlanner(simdir)
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
            error("Infeasible solution for simdir $i: $(simdirs[i])")
        end
        if termination_status(jump_model) != MOI.OPTIMAL
            println("Warning: Termination status not optimal for simdir $i")
        end
        
        export_results!(planner)
        jump_models[i] = planner.model.jump_model
        datas[i] = planner.data
        objectives[i] = objective_value(jump_model)
    end
    
    # Compute weighted average (order is preserved from create_simdirs)
    second_lb_val = sum(objectives[i] * rep_prob[i] for i in 1:n_sims)
    
    # Write to CSV
    df = DataFrame(second_stage_lb = second_lb_val)
    CSV.write(joinpath(second_lb_dir, "second_stage_lb.csv"), df)
    
    return second_lb_val, jump_models, datas
end

function compute_benders_lb(superdir::String; force_lb::Bool=false)
    """
    Computes a valid global lower bound by summing first-stage and second-stage LBs
    """
    lb_dir = mkpath(joinpath(superdir, "lower_bound"))
    first_lb_dir = mkpath(joinpath(lb_dir, "first_stage"))
    second_lb_dir = mkpath(joinpath(lb_dir, "second_stage"))
    first_lb_file = joinpath(first_lb_dir, "first_stage_lb.csv")
    second_lb_file = joinpath(second_lb_dir, "second_stage_lb.csv")
    
    # Compute or load first-stage LB
    if force_lb || !isfile(first_lb_file)
        first_stage_lb, = compute_first_stage_lb(superdir, first_lb_dir)
    else
        first_stage_lb = CSV.read(first_lb_file, DataFrame)[1, :first_stage_lb]
    end
    # Compute or load second-stage LB
    if force_lb || !isfile(second_lb_file)
        second_stage_lb, = compute_second_stage_lb(superdir, second_lb_dir)
    else
        second_stage_lb = CSV.read(second_lb_file, DataFrame)[1, :second_stage_lb]
    end
    
    benders_lb = first_stage_lb + second_stage_lb
    println("Benders Lower Bound: $benders_lb (first: $first_stage_lb, second: $second_stage_lb)")
    
    return benders_lb
end