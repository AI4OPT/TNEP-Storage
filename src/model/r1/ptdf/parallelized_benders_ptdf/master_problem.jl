using JuMP
using Gurobi
using Serialization

function save_master_problem(master, simdir; filename="master_problem.mps")
    """
    Save the master problem model to file for warm-starting
    """
    master_path = joinpath(simdir, "master_models")
    mkpath(master_path)
    
    # Save the MPS file
    mps_file = joinpath(master_path, filename)
    write_to_file(master, mps_file)
    
    # Save extension data (iterations tracking, etc.)
    extension_file = joinpath(master_path, replace(filename, ".mps" => "_extensions.jls"))
    extension_data = Dict(
        "y_raw" => master.ext[:y_raw],
        "y_eval" => master.ext[:y_eval],
        "y_core" => master.ext[:y_core],
        "iter" => master.ext[:iter],
        "gap" => master.ext[:gap],
        "stabilization_lambda" => master.ext[:stabilization_lambda],
        "stabilization_lambda_decr" => master.ext[:stabilization_lambda_decr],
        "date_weights" => master.ext[:date_weights]
    )
    
    serialize(extension_file, extension_data)
    
    println("[DEBUG] Master problem saved to $mps_file")
    println("[DEBUG] Extensions saved to $extension_file")
end

function load_or_make_master_problem(simdir, data, date_weights; filename="master_problem.mps")
    """
    Load a previously saved master problem for warm-starting
    """
    master_path = joinpath(simdir, "master_models")
    mps_file = joinpath(master_path, filename)
    extension_file = joinpath(master_path, replace(filename, ".mps" => "_extensions.jls"))
    
    if !isfile(mps_file)
        println("[DEBUG] No saved master problem found at $mps_file")
        master, y, theta = define_master_ptdf(data, date_weights)
        return master, y, theta
    end
    
    # Read the model from file (this creates a new model)
    optimizer = Gurobi.Optimizer
    master = read_from_file(mps_file)
    set_optimizer(master, optimizer)
    
    if !get(data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])
    end
    
    # Restore extension data
    if isfile(extension_file)
        extension_data = deserialize(extension_file)
        master.ext[:y_raw] = extension_data["y_raw"]
        master.ext[:y_eval] = extension_data["y_eval"]
        master.ext[:y_core] = extension_data["y_core"]
        master.ext[:iter] = extension_data["iter"]
        master.ext[:gap] = extension_data["gap"]
        master.ext[:stabilization_lambda] = extension_data["stabilization_lambda"]
        master.ext[:stabilization_lambda_decr] = extension_data["stabilization_lambda_decr"]
        master.ext[:date_weights] = extension_data["date_weights"]
        println("[DEBUG] Restored extension data with $(length(extension_data["y_raw"])) previous iterations")
    else
        error("Loaded master JuMP model but missing master extension files")
    end
    
    # Reconstruct variable references based on the loaded model
    # This assumes the variable ordering is preserved
    E = length(data["branch"])
    N = length(data["bus"])
    
    # Get variable references from the loaded model
    all_vars = all_variables(master)
    
    # Reconstruct gamma variables (first E variables)
    gamma_vars = all_vars[1:E]
    master[:gamma] = gamma_vars
    
    # Reconstruct s_power variables (next N variables)
    s_power_vars = all_vars[E+1:E+N]
    master[:s_power] = s_power_vars
    
    # Reconstruct s_energy variables (next N variables)
    s_energy_vars = all_vars[E+N+1:E+2*N]
    master[:s_energy] = s_energy_vars
    
    # Reconstruct theta variable (last variable)
    theta_vars = all_vars[E+2*N+1:end]
    master[:theta] = theta_vars
    
    y = (master[:gamma], master[:s_power], master[:s_energy])
    
    println("[DEBUG] Successfully loaded master problem from $mps_file")
    return master, y, master[:theta]
end

function define_master_ptdf(data::Dict{String, Any}, date_weights::Dict{Int, Tuple{String, Float64}})
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

    # Create extension for tracking iterations of y_raw, y_eval, and y_core
    master.ext[:y_raw] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_eval] = Vector{Vector{Vector{Float64}}}()
    master.ext[:y_core] = Vector{Vector{Vector{Float64}}}()

    # Create extensions for dynamic lambda adjustments
    master.ext[:gap] = Vector{Float64}()
    master.ext[:stabilization_lambda] = Vector{Float64}()
    master.ext[:stabilization_lambda_decr] = Vector{Float64}()
    master.ext[:iter] = 0

    # Create extension for date_weights
    master.ext[:date_weights] = date_weights

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
    if haskey(data["param"], "relaxed_first_stage") && data["param"]["relaxed_first_stage"] == true
        JuMP.@variable(master, 0 <= gamma[a=1:E] <= K)
    else
        JuMP.@variable(master, 0 <= gamma[a=1:E] <= K, Int)
    end

    # binary variable for installation of storage
    # JuMP.@variable(master, sigma[i=1:N], Bin)

    # power rating of storage
    JuMP.@variable(master, s_power[i=1:N] >= 0)

    # energy rating of storage
    JuMP.@variable(master, s_energy[i=1:N] >= 0)

    # subproblem objective(s)
    JuMP.@variable(master, theta[r=1:R] >= 0)

    #
    #   II. Constraints
    #
    candidate_branches = get_candidate_branches(data)
    
    # Check if previous investments exist
    if haskey(data["param"], "previous_investment_dir")
        prev_dir = data["param"]["previous_investment_dir"]
        add_prev_upgrades(master, data, gamma, s_energy, prev_dir, candidate_branches)
    end

    # Separately, force non-candidates to zero
    if candidate_branches !== nothing
        @constraint(master, 
            force_zero[a in rate_a_nonzero; !(a in candidate_branches)],
            gamma[a] == 0
        )
    end

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(master, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # max energy rating
    JuMP.@constraint(master, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= data["param"]["max_energy_rating"]
    )

    # max power rating
    JuMP.@constraint(master, 
        installed_power_ub[i in 1:N],
        s_power[i] <= data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    JuMP.@constraint(master, 
        short_duration[i in 1:N],
        s_energy[i] == 4.0 * s_power[i]
    )
    #
    #   III. Objective
    #
    # Is set in the Bender's loop
    y = gamma, s_power, s_energy

    return master, y, theta
end

function update_master_objective!(superdir, master, data, y, theta)
    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    date_weights = master.ext[:date_weights]

    #
    #   III. Objective
    #
    gamma, s_power, s_energy = y

    # Start with the base objective terms
    obj_expr = (sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            sum(theta[r] * date_weights[r][2] for r in 1:R))

    # L2 (quadratic) regularization
    if haskey(data["param"], "reg_penalty")
        gamma_reg, s_power_reg, s_energy_reg = compute_regularization_point(superdir, master, data)
        obj_expr += data["param"]["reg_penalty"] * sum((s_power[i] - s_power_reg[i])^2 for i in 1:N)
        if haskey(data["param"], "trans_reg_penalty")
            obj_expr += data["param"]["trans_reg_penalty"] * sum((gamma[a] - gamma_reg[a])^2 for a in 1:E)
        end
    end

    # Update the objective    
    JuMP.@objective(master, Min, obj_expr)
end

function add_benders_cut_ptdf(master, theta, duals, y, y_val, phi_val)
    # Unpack y variables and duals
    gamma, s_power, s_energy = y
    dual_gamma, dual_power = duals
    gamma_val, s_power_val, s_energy_val = y_val

    # Add Benders cut
    @constraint(master, 
        theta >= phi_val + 
        sum(dual_gamma[a] * (gamma[a] - gamma_val[a]) for a=1:length(gamma)) +
        sum(dual_power[i] * (s_power[i] - s_power_val[i]) for i=1:length(s_power))
    )
end

function add_prev_upgrades(master, data, gamma, s_energy, prev_dir, candidate_branches=nothing; tolerance=1e-5)
    E = length(data["branch"])
    N = length(data["bus"])

    trans_file = joinpath(prev_dir, "line_investments.csv")
    storage_file = joinpath(prev_dir, "storage_investments.csv")
    trans_df = CSV.read(trans_file, DataFrame)
    storage_df = CSV.read(storage_file, DataFrame)

    if candidate_branches !== nothing
        valid_branches = [a for a in 1:E if (a in candidate_branches) && (trans_df[a, :Upgrade_Lvl] > tolerance)]
    else
        # Original behavior if no candidate set provided
        valid_branches = [a for a in 1:E if trans_df[a, :Upgrade_Lvl] > tolerance]
    end

    @constraint(master,
        old_gamma[a in valid_branches],
        gamma[a] >= trans_df[a, :Upgrade_Lvl]
    )
    @constraint(master,
        old_s_energy[i in 1:N, storage_df[i, :Storage_Energy] > tolerance],
        s_energy[i] >= storage_df[i, :Storage_Energy]
    )
end

function get_candidate_branches(data, tolerance=1e-5)
    if !get(data["param"], "line_candidates", false)
        return nothing
    end
    
    dates = data["param"]["dates"] # Read data params
    initial_optima_dir = data["param"]["initial_optima_dir"]
    
    branch_rep_upgrades = Dict{Int, Vector{Float64}}() # Initialize storage for all branch data
    
    for branch_idx in 1:length(data["branch"]) # Initialize with empty vectors for each branch
        branch_rep_upgrades[branch_idx] = Float64[]
    end
    
    for a_date in dates # Collect rep day data
        rep_file = joinpath(initial_optima_dir, a_date, "output", "line_investments.csv")
        rep = CSV.read(rep_file, DataFrame)
        
        for row in eachrow(rep) # Add each branch's upgrade level to its collection
            push!(branch_rep_upgrades[row.Branch_Index], row.Upgrade_Lvl)
        end
    end

    filtered_dict = Dict(k => v for (k, v) in branch_rep_upgrades if any(abs(val) > tolerance for val in v))
    return collect(keys(filtered_dict))
end



