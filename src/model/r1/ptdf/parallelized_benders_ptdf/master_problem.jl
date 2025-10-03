using JuMP
using Gurobi
using Serialization

function define_master_ptdf(superdir, data::Dict{String, Any}, date_weights::Dict{Int, Tuple{String, Float64}})
    # Initialize model
    optimizer = Gurobi.Optimizer
    master = JuMP.Model(optimizer)
    
    if !get(data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])
    end

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Create extension for tracking iterations of y_trust
    master.ext[:y_trust] = Vector{Vector{Vector{Float64}}}()
    master.ext[:trust_transmission] = Vector{Float64}()
    master.ext[:trust_sizing] = Vector{Float64}()
    master.ext[:trust_siting] = Vector{Float64}()

    # Create extensions for dynamic lambda adjustments
    master.ext[:gap] = Vector{Float64}()
    master.ext[:stabilization_lambda] = Vector{Float64}()
    master.ext[:stabilization_lambda_decr] = Vector{Float64}()
    master.ext[:iter] = 0
    master.ext[:data] = data
    master.ext[:total_ue] = [Inf]
    master.ext[:total_obj] = [Inf]
    master.ext[:last_y_val] = [zeros(E), zeros(N)]

    # Create extension for date_weights, superdir
    master.ext[:date_weights] = date_weights
    master.ext[:superdir] = superdir

    # Get sets of branches with/without thermal limits
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)

    # Pre-compute useful mappings that don't change during the solution process
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))
    
    # Branch metadata
    from_bus = Dict(parse(Int, a) => data["branch"]["$a"]["f_bus"] for a in keys(data["branch"]))
    to_bus = Dict(parse(Int, a) => data["branch"]["$a"]["t_bus"] for a in keys(data["branch"]))

    #
    #   I. Variables
    #

    # investment level of capacity upgrade
    JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)

    # energy rating of storage
    JuMP.@variable(master, s_energy_int[i=1:N] >= 0, Int)

    # subproblem objective(s)
    JuMP.@variable(master, theta[r=1:R] >= 0)

    # subproblem load shed
    JuMP.@variable(master, ue_sum >= 0)

    master.ext[:warmstart] = get(data["param"], "warmstart", false)
    master.ext[:stabilization] = get(data["param"], "stabilization", false)
    master.ext[:level_set] = get(data["param"], "level_set", false)

    if master.ext[:stabilization] == "trust_region"
        JuMP.@variable(master, abs_diff[1:N] >= 0)
        JuMP.@variable(master, trans_abs_diff[1:E] >= 0)
    end

    #
    #   II. Constraints
    #
    if master.ext[:warmstart]
        JuMP.@constraint(master, 
            no_load_shed,
            ue_sum <= 0
        )
    end    

    # Check if previous investments exist
    if haskey(data["param"], "previous_investment_dir")
        prev_dir = data["param"]["previous_investment_dir"]
        add_prev_upgrades(master, data, prev_dir)
    end

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(master, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # energy rating maximum
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy_int[i] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
    )

    master.ext[:over_invested_point] = compute_superset_core_point(superdir)

    # fix optimistic over-invested transmission
    gamma_over_invested = master.ext[:over_invested_point][1]
    if get(data["param"], "over_invested_transmission", false)
        @constraint(master,
            over_invest_gamma[a in 1:E],
            gamma[a] >= gamma_over_invested[a]
        )
    end

    # storage candidates
    if get(data["param"], "storage_cand", false)
        storage_upgrades = master.ext[:over_invested_point][2]
        storage_non_indices = findall(x -> x == 0, storage_upgrades)

        JuMP.@constraint(master, 
            energy_cand[i in storage_non_indices],
            s_energy_int[i] == 0
        )
    end

    # line candidates
    if get(data["param"], "line_cand", false)
        line_upgrades = master.ext[:over_invested_point][1]
        line_non_indices = findall(x -> x == 0, line_upgrades)

        JuMP.@constraint(master, 
            line_cand[a in line_non_indices],
            gamma[a] == 0
        )
    end

    #
    #   III. Objective
    #
    # Start with the base objective terms
    obj_expr = sum(s_energy_int[i] * data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            sum(theta[r] * date_weights[r][2] for r in 1:R) + ue_sum * data["param"]["operational_weight"] * data["param"]["under_served_penalty"]

    y = gamma, s_energy_int

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)
    master.ext[:obj_expr] = obj_expr

    return master, y, theta
end

function add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
    # Unpack y variables and duals
    gamma, s_energy_int = y
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_int_val = y_val

    # Add Benders cut
    @constraint(master, 
        theta >= phi_val + 
        sum(dual_gamma[a] * (gamma[a] - gamma_val[a]) for a=1:length(gamma)) +
        sum(dual_energy[i] * (s_energy_int[i] - s_energy_int_val[i]) for i=1:length(s_energy_int))
    )
end

function add_prev_upgrades(master, data, prev_dir)
    E = length(data["branch"])
    N = length(data["bus"])
    gamma = master[:gamma]
    s_energy_int = master[:s_energy_int]

    trans_file = joinpath(prev_dir, "line_investments.csv")
    storage_file = joinpath(prev_dir, "storage_investments.csv")
    trans_df = CSV.read(trans_file, DataFrame)
    storage_df = CSV.read(storage_file, DataFrame)

    nonzero_trans_indices = findall(x -> x > 0, trans_df[:, :Upgrade_Lvl])
    nonzero_storage_indices = findall(x -> x > 0, storage_df[:, :Storage_Energy])

    @constraint(master,
        old_gamma[a in nonzero_trans_indices],
        gamma[a] >= trans_df[a, :Upgrade_Lvl]
    )
    @constraint(master,
        old_s_energy[i in nonzero_storage_indices],
        s_energy_int[i] >= storage_df[i, :Storage_Energy]
    )
end

function add_trust_region!(master)
    if master.ext[:stabilization] != "trust_region"
        return
    end

    data = master.ext[:data]
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    gamma, s_energy_int = master[:gamma], master[:s_energy_int]
    abs_diff, trans_abs_diff = master[:abs_diff], master[:trans_abs_diff]
    
    if master.ext[:iter] == 0
        # warm start initially to the exact over invested point
        gamma_core, s_energy_int_core = zeros(E), zeros(N)
        if master.ext[:warmstart]
            gamma_core, s_energy_int_core = master.ext[:over_invested_point]
        end
        y_core = [gamma_core, s_energy_int_core]
        push!(master.ext[:y_trust], y_core)
        master.ext[:l1_radius] = [1]
        
        y_trust = master.ext[:y_trust][end]
        gamma_trust, s_energy_int_trust = y_trust

        @constraint(master,
            trans_abs_diff_ub[a in 1:E],
            trans_abs_diff[a] >= gamma[a] - gamma_trust[a])
        @constraint(master,
            trans_abs_diff_lb[a in 1:E],
            trans_abs_diff[a] >= gamma_trust[a] - gamma[a])

        @constraint(master,
            abs_diff_ub[i in 1:N],
            abs_diff[i] >= s_energy_int[i] - s_energy_int_trust[i])
        @constraint(master,
            abs_diff_lb[i in 1:N],
            abs_diff[i] >= s_energy_int_trust[i] - s_energy_int[i])

        @constraint(master,
            trans_abs_diff_total,
            sum(trans_abs_diff[a] for a in 1:E) == 0)
        @constraint(master,
            abs_diff_total,
            sum(abs_diff[i] for i in 1:N) == 0)
        return
    else
        remove_trust_region!(master)
    end
    
    # Move these assignments outside the if-else block to ensure they're always defined
    y_trust = master.ext[:y_trust][end]
    gamma_trust, s_energy_int_trust = y_trust

    @constraint(master,
        trans_abs_diff_ub[a in 1:E],
        trans_abs_diff[a] >= gamma[a] - gamma_trust[a])
    @constraint(master,
        trans_abs_diff_lb[a in 1:E],
        trans_abs_diff[a] >= gamma_trust[a] - gamma[a])
    
    @constraint(master,
        abs_diff_ub[i in 1:N],
        abs_diff[i] >= s_energy_int[i] - s_energy_int_trust[i])
    @constraint(master,
        abs_diff_lb[i in 1:N],
        abs_diff[i] >= s_energy_int_trust[i] - s_energy_int[i])

    @constraint(master,
            trans_abs_diff_total,
            sum(trans_abs_diff[a] for a in 1:E) <= master.ext[:l1_radius][end])
    @constraint(master,
            abs_diff_total,
            sum(abs_diff[i] for i in 1:N) == 2)
end

function remove_trust_region!(master)
    stab_method = master.ext[:stabilization]

    if stab_method == "trust_region"
        constraint_names = [
            :trans_abs_diff_ub, :trans_abs_diff_lb, :abs_diff_ub, :abs_diff_lb, :abs_diff_total, :trans_abs_diff_total
        ]
    elseif stab_method == "boxstep"
        constraint_names = [
        :trans_box_ub, :trans_box_lb, :stor_box_lb, :stor_box_ub
        ]
    else
        return
    end

    obj_dict = object_dictionary(master)

    for name in constraint_names
        if haskey(obj_dict, name)
            try
                constraint_obj = obj_dict[name]

                # JuMP.delete works on both single constraints and arrays
                JuMP.delete(master, constraint_obj)
                JuMP.unregister(master, name)
                println("Removed constraint: $name")

            catch e
                println("Warning: Could not remove $name: $e")
                # Force unregister
                try
                    JuMP.unregister(master, name)
                catch
                end
            end
        else
            println("Constraint $name not found in model")
        end
    end
end

function add_level_set!(master, current_obj)
    if master.ext[:level_set] != true
        return
    end
    remove_level_set!(master)

    obj_expr = master.ext[:obj_expr]
    @constraint(master,
        level_set,
        obj_expr <= current_obj
    )
end

function remove_level_set!(master)
    obj_dict = object_dictionary(master)
    name = :level_set
    if haskey(obj_dict, name)
        try
            constraint_obj = obj_dict[name]

            # JuMP.delete works on both single constraints and arrays
            JuMP.delete(master, constraint_obj)
            JuMP.unregister(master, name)
            println("Removed constraint: $name")

        catch e
            println("Warning: Could not remove $name: $e")
            # Force unregister
            try
                JuMP.unregister(master, name)
            catch
            end
        end
    else
        println("Constraint $name not found in model")
    end
end
