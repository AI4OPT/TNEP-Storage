using Gurobi
using JSON

include("../../helpers/setup_simdir.jl")
include("add_params_profiles.jl")
include("export_model.jl")
include("create_model.jl")

function run_model(simdir)
    setup_simdir(simdir)

    # load data and add params and profiles
    data = add_params_profiles(simdir)

    # create the model
    optimizer = Gurobi.Optimizer
    model = create_model(data, optimizer)

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