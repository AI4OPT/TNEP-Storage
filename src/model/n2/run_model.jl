using Gurobi
using JSON

include("../../helpers/setup_simdir.jl")
include("create_model.jl")
include("add_params_profiles.jl")
include("../n1/export_model.jl")
include("decarbonization.jl")

function run_model(simdir; prev_model=nothing) # set add_ts to be true if anything other than the params in the powermodels file needs to be changed
    setup_simdir(simdir)

    data = add_params_profiles(simdir)
    data = update_decarbonization(simdir, data)
    json_data = JSON.json(data)
    open(joinpath(simdir, "data.json"), "w") do file
        JSON.print(file, data)
    end

    # create the model
    optimizer = Gurobi.Optimizer
    model = create_model_n2(data, optimizer, prev_model=prev_model)

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
    export_flow(simdir, model, data)
end

function run_sequential(seqsimdir, ending_year)
    config_file = joinpath(seqsimdir, "config.toml")
    toml_data = TOML.parsefile(config_file)
    current_year = toml_data["decarbonization_year"]

    # make the starting directory
    simdir = joinpath(seqsimdir, "$current_year")
    if !isdir(simdir)
        mkdir(simdir)
    end

    # copy the toml to the new directory
    open(joinpath(simdir, "config.toml"), "w") do file
        TOML.print(file, toml_data)
    end

    current_model, current_data = run_model(simdir)

    while current_year != ending_year
        config_file = joinpath(simdir, "config.toml")
        toml_data = TOML.parsefile(config_file)
        toml_data["decarbonization_year"] += 1 # TODO CHECK IF THIS IS INT

        # update current_year and simdir
        current_year = toml_data["decarbonization_year"]
        simdir = joinpath(seqsimdir, "$current_year")
        if !isdir(simdir)
            mkdir(simdir)
        end

        # copy the toml to the new directory
        open(joinpath(simdir, "config.toml"), "w") do file
            TOML.print(file, toml_data)
        end

        current_model, current_model = run_model(simdir, prev_model=current_model)
    end
end

