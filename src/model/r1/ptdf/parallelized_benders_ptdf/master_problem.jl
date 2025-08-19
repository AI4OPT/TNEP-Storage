using JuMP
using Gurobi
using Serialization

function define_master_ptdf(superdir, data::Dict{String, Any}, date_weights::Dict{Int, Tuple{String, Float64}};
    max_ptdf_iterations::Int=256,
    max_ptdf_per_iteration::Int=32,
    ptdf_tol::Float64=1e-6)

    embed_date = data["param"]["embed_date"]
    
    # Initialize model
    optimizer = Gurobi.Optimizer
    simdir = joinpath(superdir, embed_date)
    sub_data = JSON.parsefile(joinpath(simdir, "data.json"))
    master = setup_ptdf_model(simdir, sub_data, optimizer, max_ptdf_iterations, max_ptdf_per_iteration, ptdf_tol)
    
    if !get(data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])
    end
    
    warm_start_ptdf_constraints!(simdir, master, sub_data)
    
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
    
    # Initialize sets - use data (master data) for dimensions, not sub_data
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    
    # Add theta variables for second stage costs (one for each representative period)
    @variable(master, theta[r in 1:R])
    
    # Get first-stage variables from the master model (these were created in setup_ptdf_model)
    gamma = master[:gamma]
    s_power = master[:s_power]
    s_energy = master[:s_energy]
    
    # Get operational variables for the embedded subproblem (these exist for the embedded date)
    pg = master[:pg]
    ue = master[:ue]
    ch = master[:ch]
    dis = master[:dis]
    
    #
    # III. Objective
    #
    # First part of objective is the investments
    first_stage_obj_expr = (
        sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
        sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E)
    )
    
    operational_weight = get(data["param"], "operational_weight", 1)
    
    # Embedded subproblem objective (for the specific embed_date)
    # Find which representative period corresponds to embed_date
    embed_rep = findfirst(r -> date_weights[r][1] == embed_date, 1:R)
    if embed_rep === nothing
        error("embed_date $embed_date not found in date_weights")
    end
    
    embed_obj = (
        sum(
            sum(compute_gen_cost(pg[1, g, t], data["gen"]["$g"]) for g=1:G) +
            sum(data["param"]["under_served_penalty"] * ue[1, i, t] for i=1:N) +
            sum(get(data["param"], "storage_operation_cost", 0.0) * (ch[1,i,t] + dis[1,i,t]) for i=1:N)
            for t=1:T
        ) * operational_weight
    )
    
    # Constraint linking theta variable to embedded objective
    @constraint(master, theta[embed_rep] >= embed_obj)
    
    # Second stage objective (sum over all representative periods)
    second_stage_obj_expr = sum(theta[r] * date_weights[r][2] for r in 1:R)
    
    # Total objective
    obj_expr = first_stage_obj_expr + second_stage_obj_expr
    @objective(master, Min, obj_expr)
    
    # Return first-stage variables as a tuple for easy access
    y = (gamma, s_power, s_energy)
    
    return master, y, theta
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



