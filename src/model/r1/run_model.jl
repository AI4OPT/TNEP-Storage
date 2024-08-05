using Gurobi
using JSON

include("../../helpers/convert_units.jl")
include("../../helpers/setup_simdir.jl")
include("add_params_profiles.jl")
include("../n1/export_model.jl")
include("create_model.jl")
include("../n2/decarbonization.jl")
include("create_summary.jl")

function run_model(simdir; prev_simdir=nothing, timeout=84600)
    setup_simdir(simdir)

    data = add_params_profiles(simdir)
    data = update_decarbonization(simdir, data)
    data = convert_units(data)
    json_data = JSON.json(data)
    open(joinpath(simdir, "data.json"), "w") do file
        JSON.print(file, data)
    end

    # create the model
    optimizer = Gurobi.Optimizer
    model = create_model_r1(data, optimizer, prev_simdir=prev_simdir)

    # set Gurobi log location
    set_optimizer_attribute(model, "LogFile", joinpath(simdir, "gurobi_logfile.log"))
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])
    set_optimizer_attribute(model, "TimeLimit", timeout)
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        println("Termination status not optimal")
        # Export the solution even if it's suboptimal
        export_model(simdir, model, data)
        write_summary_to_csv(simdir, model, data)
        return model, data
    end

    if primal_status(model) != MOI.FEASIBLE_POINT
        error("Primal status not feasible point")
    end

    export_model(simdir, model, data)
    write_summary_to_csv(simdir, model, data)

    return model, data

end

function export_model(simdir, model, data)
    export_investments_csv(simdir, model, data)
    export_energy_csv(simdir, model, data)
    export_flow(simdir, model, data)
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
