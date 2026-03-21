function compute_second_stage_costs(superdir, data, benders_iter)
    dates = data["param"]["dates"]
    date_cores = [x[6:end] for x in dates]
    probs = data["param"]["representative_prob"]
    years = data["param"]["years"]

    date_weights = Dict()
    for i in 1:length(date_cores)
        date_weights[date_cores[i]] = probs[i]
    end

    for year in years
        year_opex = 0
        for date in date_cores
            simdir = joinpath(superdir, string(year) * "-" * date)
            df = CSV.read(joinpath(simdir, "output", "benders_progress.csv"), DataFrame)
            phi_val = df[benders_iter, :phi_val]
            year_opex += (phi_val * date_weights[date])
        end
        println("Year $year: OpEx = $year_opex")
    end
end

function compute_first_stage_costs(superdir, data, benders_iter)
    years = data["param"]["years"]

    for year in years
        line_df = CSV.read(joinpath(superdir, "benders_output", "line_investments_" * string(benders_iter) * "_" * string(year) * ".csv"), DataFrame)
        stor_df = CSV.read(joinpath(superdir, "benders_output", "storage_investments_" * string(benders_iter) * "_" * string(year) * ".csv"), DataFrame)

        gamma_val = line_df[:, :Upgrade_Lvl]
        sigma_val = stor_df[:, :Storage_Energy]

        stor_cost = sum(sigma_val) * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"]
        line_cost = sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * 
                    data["branch"]["$a"]["distance"] * gamma_val[a] for a in 1:length(gamma_val))

        println("$year, first-stage cost: $(stor_cost + line_cost)")
        println("$year, trans cost: $(line_cost)")
        println("$year, stor cost: $(stor_cost)")
    end
end

function compute_first_stage_costs(dir, data)
    line_df = CSV.read(joinpath(dir, "output", "line_investments.csv"), DataFrame)
    stor_df = CSV.read(joinpath(dir, "output", "storage_investments.csv"), DataFrame)

    gamma_val = line_df[:, :Upgrade_Lvl]
    sigma_val = stor_df[:, :Storage_Energy]

    stor_cost = sum(sigma_val) * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"]
    line_cost = sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * 
                data["branch"]["$a"]["distance"] * gamma_val[a] for a in 1:length(gamma_val))

    println("$year, first-stage cost: $(stor_cost + line_cost)")
    println("$year, trans cost: $(line_cost)")
    println("$year, stor cost: $(stor_cost)")
end
