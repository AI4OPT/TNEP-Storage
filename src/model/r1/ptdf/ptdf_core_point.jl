function compute_core_point(data::Dict{String, Any}; transmission_eps=150, storage_eps=1.5)
    # Create a copy of the master problem but with relaxed integer variables
    master_cp, y_cp, theta_cp = define_master_ptdf(data)
    
    # Unpack the variables
    gamma_cp, s_power_cp, s_energy_cp = y_cp
    
    # Relax integer constraints on gamma
    unset_integer.(master_cp, gamma_cp)
    
    # Add the minimum requirements for transmission and storage
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)
    
    # For lines that can be upgraded, add a small minimum total capacity
    @constraint(master_cp, 
        sum(gamma_cp[a] for a in 1:length(gamma_cp) if !(a in rate_a_zero)) >= transmission_eps
    )
    
    # If storage_eps is provided, add constraint for minimum storage
    if storage_eps > 0
        @constraint(master_cp,
            sum(s_power_cp[i] for i in 1:length(s_power_cp)) >= storage_eps
        )
    end
    
    # Initially, the theta value from the master problem formulation may not 
    # be meaningful without cuts, so give it a reasonable lower bound
    set_lower_bound(theta_cp, 0)
    
    # Solve the relaxed problem with the original objective
    optimize!(master_cp)
    
    if termination_status(master_cp) != MOI.OPTIMAL
        error("Core point computation did not solve to optimality: $(termination_status(master_cp))")
    end
    
    # Extract solution values
    gamma_val = value.(gamma_cp)
    s_power_val = value.(s_power_cp)
    s_energy_val = value.(s_energy_cp)
    
    # Round gamma to integers for use in the master problem
    gamma_val_rounded = round.(Int, gamma_val)
    
    # Ensure constraints are satisfied after rounding
    for a in 1:length(gamma_val_rounded)
        if a in rate_a_zero
            gamma_val_rounded[a] = 0
        elseif gamma_val_rounded[a] > data["param"]["num_cap_upgrades_max"]
            gamma_val_rounded[a] = data["param"]["num_cap_upgrades_max"]
        elseif gamma_val_rounded[a] < 0
            gamma_val_rounded[a] = 0
        end
    end
    
    # Print summary of the core point
    println("Core Point Summary:")
    println("  Total line upgrades: $(sum(gamma_val_rounded))")
    println("  Total storage power: $(sum(s_power_val)) MW")
    println("  Total storage energy: $(sum(s_energy_val)) MWh")
    println("  Number of buses with storage: $(count(x -> x > 1e-4, s_power_val))")
    
    # Return the core point
    return (gamma_val_rounded, s_power_val, s_energy_val)
end

function get_over_invested_point(simdir, data)
    config_file = joinpath(simdir, "config.toml")
    config = TOML.parsefile(config_file)
    rep_date = config["dates"][1]

    if isfile(joinpath(simdir, "line_investments.csv")) && isfile(joinpath(simdir, "storage_investments.csv"))
        line_file = joinpath(simdir, "line_investments.csv")
        storage_file = joinpath(simdir, "storage_investments.csv")
    else
        output_dir = joinpath(simdir, "output")
        line_file = joinpath(output_dir, "line_investments.csv")
        storage_file = joinpath(output_dir, "storage_investments.csv")
    end

    line_inv = CSV.read(line_file, DataFrame)
    storage_inv = CSV.read(storage_file, DataFrame)

    gamma_val = line_inv[:, :Upgrade_Lvl]
    s_energy_int_val = storage_inv[:, :Storage_Energy] * 1 / data["param"]["storage_energy_size"]


    return ceil.(Int, gamma_val), ceil.(Int, s_energy_int_val)
end

function get_rep_day_core_point(simdir)
    config_file = joinpath(simdir, "config.toml")
    config = TOML.parsefile(config_file)
    rep_date = config["dates"][1]

    if isfile(joinpath(simdir, "line_investments.csv")) && isfile(joinpath(simdir, "storage_investments.csv"))
        line_file = joinpath(simdir, "line_investments.csv")
        storage_file = joinpath(simdir, "storage_investments.csv")
    else
        output_dir = joinpath(simdir, "output")
        line_file = joinpath(output_dir, "line_investments.csv")
        storage_file = joinpath(output_dir, "storage_investments.csv")
    end

    line_inv = CSV.read(line_file, DataFrame)
    storage_inv = CSV.read(storage_file, DataFrame)

    gamma_val = line_inv[:, :Upgrade_Lvl]
    s_energy_val = storage_inv[:, :Storage_Energy]

    # You can round gamma if it's binary or discrete
    # gamma_val = round.(Int, gamma_val)

    # Return the core point
    return gamma_val, s_energy_val
end
