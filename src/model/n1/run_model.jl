using Gurobi
using JSON

include("create_model.jl")
include("add_params_profiles.jl")
include("export_model.jl")

function setup_simdir(simdir)
    if !isdir(joinpath(simdir, "output"))
        mkdir(joinpath(simdir, "output"))
    end

    if !isdir(joinpath(simdir, "visual"))
        mkdir(joinpath(simdir, "visual"))
    end
end

function run_model(simdir, add_ts::Bool=false) # set add_ts to be true if anything other than the params in the powermodels file needs to be changed
    setup_simdir(simdir)

    if add_ts
        add_params_profiles(simdir)
    end

    p_data_file = joinpath(simdir, "ps_profiled_data.json")
    data = JSON.parsefile(p_data_file)
    # re-add params, in case any parameters have been changed
    data = add_params(simdir, data)

    # create the model
    optimizer = Gurobi.Optimizer
    model = create_model_n1(data, optimizer)

    # set Gurobi log location
    set_optimizer_attribute(model, "LogFile", joinpath(simdir, "gurobi_logfile.log"))
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])
    optimize!(model)

    if termination_status(model) != OPTIMAL 
        println("Termination status not optimal")
        return
    end
            
    if primal_status(model) != FEASIBLE_POINT
        error("Primal status not feasible point")
    end

    export_model(simdir, model, data)

    return model, data
end

function export_model(simdir, model, data)
    export_investments_csv(simdir, model, data)
    export_energy_csv(simdir, model, data)
end