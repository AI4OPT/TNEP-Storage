function get_rep_day_core_point(simdir)
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
    s_power_val = storage_inv[:, :Storage_Power]
    s_energy_val = storage_inv[:, :Storage_Energy]

    # You can round gamma if it's binary or discrete
    # gamma_val = round.(Int, gamma_val)

    # Return the core point
    return gamma_val, s_power_val, s_energy_val
end

function compute_regularization_point(superdir, master, data)
    iter = master.ext[:iter]

    # Compute the reg point
    gamma_reg, s_power_reg, s_energy_reg = get_rep_day_core_point(superdir)

    if iter > 0
        # Get the old eval point
        y_eval_old = master.ext[:y_eval][iter]
        gamma_reg, s_power_reg, s_energy_reg = y_eval_old
    end

    return gamma_reg, s_power_reg, s_energy_reg
end

function get_stabilization_shift(master, data)
    iter = master.ext[:iter]

    # Get lambda parameters
    stabilization_lambda = get(data["param"], "stabilization_lambda", [0.0, 0.0, 0.0])
    starter_lambda = stabilization_lambda[1]
    decrement = stabilization_lambda[2]
    min_lambda = stabilization_lambda[3]
    lambda = nothing
    decr = nothing
    # Get core shift phi parameters
    core_shift_phi = get(data["param"], "core_shift_phi", [0.0, 0.0])
    main_phi = core_shift_phi[1]
    low_phi = core_shift_phi[2]
    phi = nothing

    if isempty(master.ext[:gap]) # i.e. iter == 0
        push!(master.ext[:stabilization_lambda], starter_lambda)
        push!(master.ext[:stabilization_lambda_decr], decrement)
        return starter_lambda, low_phi
    else
        if haskey(data["param"], "dynamic_lambda")
            # Handle dynamic lambda
            gap = master.ext[:gap][end]
            threshold, EXTREME_GAP = data["param"]["dynamic_lambda"]
            old_lambda = master.ext[:stabilization_lambda][end]
            old_decr = master.ext[:stabilization_lambda_decr][end]

            if length(master.ext[:gap]) > 5 && all(master.ext[:gap][end-4:end] .> threshold * 0.5)
                lambda = old_lambda + decrement / 4
                decr = decrement / 8
                phi = low_phi 
            elseif gap > threshold
                # increment lambda back (retreat into stabilized region)
                lambda = min(1.0, old_lambda + old_decr / 2)
                # reduce the size of the decrementor
                decr = old_decr / 4
                phi = low_phi
            else
                # decrease lambda (explore new areas more)
                lambda = max(0.0, old_lambda - old_decr)
                # increase the size of the decrementor
                decr = old_decr * 2
                phi = main_phi
            end
        else
            # Handle linear lambda decrease
            lambda = max(min_lambda, starter_lambda - iter * decrement)
            phi = main_phi
        end
    end

    push!(master.ext[:stabilization_lambda], lambda)
    push!(master.ext[:stabilization_lambda_decr], decr)
    return lambda, phi
end

function compute_eval_core_points(superdir, master, data; trans_ratio=0.005, er_ratio=0.0005)
    iter = master.ext[:iter]

    trans_ratio = get(data["param"], "trans_perturb_ratio", trans_ratio)
    er_ratio = get(data["param"], "er_perturb_ratio", er_ratio)

    # Compute stabilization lambda
    stab_lambda, shift_phi = get_stabilization_shift(master, data)

    # core point perturbations
    gamma_eps = trans_ratio * data["param"]["num_cap_upgrades_max"]
    s_energy_eps = er_ratio * data["param"]["max_energy_rating"]
    s_power_eps =  s_energy_eps * 1/4

    y_core = nothing
    if iter == 0
        # Compute the core point
        gamma_core, s_power_core, s_energy_core = get_rep_day_core_point(superdir)
        gamma_core = Float64.(gamma_core)
    
        gamma_core = max.(gamma_core, gamma_eps)
        s_power_core = max.(s_power_core, s_power_eps)
        s_energy_core = 4.0 * s_power_core

        y_core = [gamma_core, s_power_core, s_energy_core]
        push!(master.ext[:y_core], y_core)

    else
        y_core_old = master.ext[:y_core][iter]
        y_eval_old = master.ext[:y_eval][iter]

        # Compute new shifted core point (convex combo of old core with old eval point)
        y_core = shift_phi * y_eval_old + (1 - shift_phi) * y_core_old
        push!(master.ext[:y_core], y_core)
    end

    # Compute the eval point (convex combo of raw with the core point)
    y_raw = master.ext[:y_raw][iter + 1]
    y_eval = stab_lambda * y_core + (1 - stab_lambda) * y_raw
    push!(master.ext[:y_eval], y_eval)

    return y_eval, y_core
end

