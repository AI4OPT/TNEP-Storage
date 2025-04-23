using Gurobi
# using CPLEX
using JSON

include("../../helpers/convert_units.jl")
include("../../helpers/setup_simdir.jl")
include("add_params_profiles.jl")
include("../n1/export_model.jl")
include("create_model.jl")
include("create_model_storage_linearized.jl")
include("create_model_simplified.jl")
include("../n2/decarbonization.jl")
include("create_summary.jl")
include("ptdf/create_model_ptdf.jl")
include("ptdf/ptdf_iterative.jl")
include("ptdf/ptdf_iterative_simplified.jl")
include("ptdf/ptdf_iterative_simplified_sorted.jl")

function run_model(simdir; timeout=84600)
    setup_simdir(simdir)

    data_path = joinpath(simdir, "data.json")
    data = nothing
    # if isfile(data_path)
    if 0 == 1
        # If file exists, read it directly
        data = JSON.parsefile(data_path)
    else
        # Otherwise, perform the other actions and write the file
        data = add_params_profiles(simdir)
        data = update_decarbonization(simdir, data)
        data = convert_units(data)
        
        # Write the JSON file once profiles are added and units fixed
        open(data_path, "w") do file
            JSON.print(file, data)
        end
    end

    # check if there are prior investments file
    line_investments = isfile(joinpath(simdir, "line_investments.csv")) ? joinpath(simdir, "line_investments.csv") : nothing
    storage_investments = isfile(joinpath(simdir, "storage_investments.csv")) ? joinpath(simdir, "storage_investments.csv") : nothing

    # create the model
    # optimizer = CPLEX.Optimizer
    optimizer = Gurobi.Optimizer

    model = nothing
    if !haskey(data["param"], "storage_linearized") || data["param"]["storage_linearized"] == false
        model = create_model_r1(data, optimizer,  
            line_investments=line_investments, 
            storage_investments=storage_investments)
    elseif data["param"]["storage_linearized"] == "ptdf"
        model = create_model_r1_ptdf_iterative(simdir, data, optimizer,  
            line_investments=line_investments, 
            storage_investments=storage_investments)
    elseif data["param"]["storage_linearized"] == "ptdf_simplified"
        model = create_model_r1_ptdf_iterative_simplified(simdir, data, optimizer,  
            line_investments=line_investments, 
            storage_investments=storage_investments)
    elseif data["param"]["storage_linearized"] == "ptdf_simplified_sorted"
        model = create_model_r1_ptdf_iterative_simplified_sorted(simdir, data, optimizer,  
            line_investments=line_investments, 
            storage_investments=storage_investments)
    else
        model = create_model_r1_sl(data, optimizer,  
            line_investments=line_investments, 
            storage_investments=storage_investments)
    end

    # set Gurobi log location
    set_optimizer_attribute(model, "LogFile", joinpath(simdir, "gurobi_logfile.log"))
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])
    set_optimizer_attribute(model, "TimeLimit", timeout)
    optimize!(model)

    if primal_status(model) != MOI.FEASIBLE_POINT
        # error("Primal status not feasible point")
        return model, data
    end

    if termination_status(model) != MOI.OPTIMAL
        println("Termination status not optimal")
        # Export the solution even if it's suboptimal
        export_model(simdir, model, data)
        write_summary_to_csv(simdir, model, data)
        return model, data
    end

    export_model(simdir, model, data)
    write_summary_to_csv(simdir, model, data)

    # save congested lines and hours
    if !haskey(data["param"], "storage_linearized") || data["param"]["storage_linearized"] == false
        save_congested(simdir, model, data)
    end

    return model, data

end

function run_sequential(seqsimdir, ending_year)
    ending_year = parse(Int, ending_year)
    config_file = joinpath(seqsimdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    current_year = toml_data["decarbonization_year"]

    while current_year <= ending_year
        # make the starting directory
        simdir = joinpath(seqsimdir, "$current_year")
        if !isdir(simdir)
            mkdir(simdir)
        end

        # create the new toml in the directory
        open(joinpath(simdir, "config.toml"), "w") do file
            TOML.print(file, toml_data)
        end

        # check if previous year directory exists, and run
        prev_year = current_year - 1
        if !isdir(joinpath(seqsimdir, "$prev_year"))
            run_model(simdir)
        else
            run_model(simdir, prev_simdir=joinpath(seqsimdir, "$prev_year"))
        end

        # increment current_year
        current_year += 1
        toml_data["decarbonization_year"] = current_year
    end
end

function display_active_constraints(model)
    for constr in all_constraints(model)
        dual_val = dual(constr)
        if abs(dual_val) > 1e-5  # Threshold for considering a constraint as active
            println("Constraint: ", constr, " is active with dual value: ", dual_val)
        end
    end
end
