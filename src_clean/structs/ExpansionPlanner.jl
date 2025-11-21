using Gurobi
using JSON

include("../../helpers/convert_units.jl")
include("../../helpers/setup_simdir.jl")
include("add_params_profiles.jl")
include("../n1/export_model.jl")
include("../n2/decarbonization.jl")
include("create_summary.jl")

function set_up_data(simdir)
    data_path = joinpath(simdir, "data.json")
    data = add_params_profiles(simdir)
    data = update_decarbonization(simdir, data)
    data = convert_units(data)
    
    # Write the JSON file once profiles are added and units fixed
    open(data_path, "w") do file
        JSON.print(file, data)
    end

    return data
end

mutable struct ExpansionPlanner
    simdir::String
    data::Dict{String, Any}
    optimizer::Type
    model::Union{OptimizationModel, Nothing}
    timeout::Int
    
    function ExpansionPlanner(simdir::String; 
                             optimizer_type=Gurobi.Optimizer,
                             timeout=84600)
        setup_simdir(simdir)
        data = set_up_data(simdir)
        new(simdir, data, optimizer_type, nothing, timeout)
    end
end

# Method to create the appropriate model based on data parameters
function create_model!(planner::ExpansionPlanner)
    data = planner.data
    optimizer = planner.optimizer
        
    if data["param"]["storage_linearized"] == "phase_angle"
        println("create_model_phase_angle")
        jump_model = create_model_phase_angle(data, optimizer)
        planner.model = PhaseAngleModel(jump_model, data, simdir)
        
    elseif data["param"]["storage_linearized"] == "ptdf_simplified_sorted"
        println("create_model_r1_ptdf_iterative_simplified_sorted_efficiency")
        jump_model = create_model_r1_ptdf_iterative_simplified_sorted_efficiency(
            planner.simdir, data, optimizer)
        planner.model = PTDFModel(jump_model, data, simdir)
    else
        error("Unknown storage_linearized value: $storage_linearized")
    end
    
    return planner.model
end

# Multiple dispatch for configure_optimizer!
function configure_optimizer!(planner::ExpansionPlanner)
    model = planner.model.jump_model
    data = planner.data
    
    set_optimizer_attribute(model, "LogFile", 
                           joinpath(planner.simdir, "gurobi_logfile.log"))
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])
    set_optimizer_attribute(model, "TimeLimit", planner.timeout)
end

# Multiple dispatch for optimize! - default behavior
function optimize!(planner::ExpansionPlanner)
    optimize!(planner.model.jump_model)
end

function export_results!(planner::ExpansionPlanner)
    model = planner.model
    export_results!(model)
end

function export_results!(model::OptimizationModel)
    # Call base export
    export_model(model.simdir, model.jump_model, model.data)
    write_summary_to_csv(model.simdir, model.jump_model, model.data)
end

# Main run function
function run_model(simdir::String; timeout=84600)
    planner = ExpansionPlanner(simdir; timeout=timeout)
    create_model!(planner)
    configure_optimizer!(planner)
    
    # This will dispatch to the right optimize! method
    optimize!(planner.model)
    
    # Check status
    jump_model = planner.model.jump_model
    if primal_status(jump_model) != MOI.FEASIBLE_POINT
        return jump_model, planner.data
    end
    
    if termination_status(jump_model) != MOI.OPTIMAL
        println("Termination status not optimal")
    end
    
    export_results!(planner)
    
    return jump_model, planner.data
end