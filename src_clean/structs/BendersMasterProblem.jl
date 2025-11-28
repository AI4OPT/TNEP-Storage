using JuMP
using Gurobi
using Serialization
using CSV, DataFrames

# Benders Master Problem struct
mutable struct BendersMasterProblem
    jump_model::JuMP.Model
    data::Dict{String, Any}
    superdir::String
    
    # Decision variables (stored for convenience)
    gamma::Vector{VariableRef}
    s_energy_int::Vector{VariableRef}
    theta::Vector{VariableRef}
    ue_sum::VariableRef
    
    # Metadata and tracking
    date_weights::Dict{Int, Tuple{String, Float64}}
    y_trust::Vector{Vector{Vector{Float64}}}
    gap::Vector{Float64}
    iter::Int
    total_ue::Vector{Float64}
    total_obj::Vector{Float64}
    last_y_val::Vector{Vector{Float64}}
    over_invested_point::Tuple{Vector{Float64}, Vector{Float64}}
    
    # Settings
    warmstart::Bool
    stabilization::Union{String, Bool}
    level_set::Bool
    
    # Dimension sizes
    R::Int  # num_representatives
    N::Int  # num buses
    E::Int  # num branches
    T::Int  # num hours
    G::Int  # num generators
    K::Int  # num_cap_upgrades_max
    
    # Inner constructor
    function BendersMasterProblem(superdir::String, data::Dict{String, Any}, 
                                 date_weights::Dict{Int, Tuple{String, Float64}})
        
        # Initialize model
        optimizer = Gurobi.Optimizer
        jump_model = JuMP.Model(optimizer)
        
        if !get(data["param"], "relaxed_first_stage", false)
            set_optimizer_attribute(jump_model, "MIPGap", data["param"]["mip_gap"])
        end
        
        # Get dimension sizes
        R = data["param"]["num_representatives"]
        N = length(data["bus"])
        E = length(data["branch"])
        T = data["param"]["num_hours"]
        G = length(data["gen"])
        K = data["param"]["num_cap_upgrades_max"]
        
        # Get rate_a sets
        rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
        rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
        
        # Settings
        warmstart = get(data["param"], "warmstart", false)
        stabilization = get(data["param"], "stabilization", false)
        level_set_flag = get(data["param"], "level_set", false)
        
        #
        #   I. Variables
        #
        
        # Investment level of capacity upgrade
        @variable(jump_model, 0 <= gamma[a=1:E] <= K, Int)
        
        # Energy rating of storage
        @variable(jump_model, s_energy_int[i=1:N] >= 0, Int)
        
        # Subproblem objectives
        @variable(jump_model, theta[r=1:R] >= 0)
        
        # Subproblem load shed
        @variable(jump_model, ue_sum >= 0)
        
        # Trust region variables if needed
        if stabilization == "trust_region"
            @variable(jump_model, abs_diff[1:N] >= 0)
            @variable(jump_model, trans_abs_diff[1:E] >= 0)
        end
        
        #
        #   II. Constraints
        #
        
        if warmstart
            @constraint(jump_model, no_load_shed, ue_sum <= 0)
        end
        
        # Check if previous investments exist
        if haskey(data["param"], "previous_investment_dir")
            prev_dir = data["param"]["previous_investment_dir"]
            add_prev_upgrades_internal!(jump_model, data, prev_dir, gamma, s_energy_int)
        end
        
        # If rate a is zero (unlimited), then don't allow upgrades
        @constraint(jump_model, 
            rate_a_zero_line_upgrade[a in rate_a_zero],
            gamma[a] == 0
        )
        
        # Energy rating maximum
        @constraint(jump_model, 
            installed_energy_ub[i in 1:N],
            s_energy_int[i] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
        )
        
        # Compute over-invested point
        over_invested_point = compute_superset_core_point(superdir)
        export_investments_csv(data, over_invested_point[1], over_invested_point[2], output_dir=joinpath(superdir, "warmstart"))
        
        # Storage candidates
        if get(data["param"], "storage_cand", false)
            storage_upgrades = over_invested_point[2]
            storage_non_indices = findall(x -> x == 0, storage_upgrades)
            
            @constraint(jump_model, 
                energy_cand[i in storage_non_indices],
                s_energy_int[i] == 0
            )
        end
        
        # Line candidates
        if get(data["param"], "line_cand", false)
            line_upgrades = over_invested_point[1]
            line_non_indices = findall(x -> x == 0, line_upgrades)
            
            @constraint(jump_model, 
                line_cand[a in line_non_indices],
                gamma[a] == 0
            )
        end
        
        #
        #   III. Objective
        #
        
        obj_expr = (
            sum(s_energy_int[i] * data["param"]["storage_energy_size"] * 
                data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * 
                data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
            sum(theta[r] * date_weights[r][2] for r in 1:R) + 
            ue_sum * data["param"]["operational_weight"] * data["param"]["under_served_penalty"]
        )
        
        @objective(jump_model, Min, obj_expr)
        
        # Store objective expression for later use
        jump_model.ext[:obj_expr] = obj_expr
        
        # Create and return the struct
        new(jump_model, data, superdir,
            gamma, s_energy_int, theta, ue_sum,
            date_weights,
            Vector{Vector{Vector{Float64}}}(),  # y_trust
            Vector{Float64}(),  # gap
            0,  # iter
            [Inf],  # total_ue
            [Inf],  # total_obj
            [zeros(E), zeros(N)],  # last_y_val
            over_invested_point,
            warmstart, stabilization, level_set_flag,
            R, N, E, T, G, K)
    end
end

# Helper function for adding previous upgrades (internal use)
function add_prev_upgrades_internal!(model::JuMP.Model, data, prev_dir, gamma, s_energy_int)
    E = length(data["branch"])
    N = length(data["bus"])
    
    trans_file = joinpath(prev_dir, "line_investments.csv")
    storage_file = joinpath(prev_dir, "storage_investments.csv")
    trans_df = CSV.read(trans_file, DataFrame)
    storage_df = CSV.read(storage_file, DataFrame)
    
    nonzero_trans_indices = findall(x -> x > 0, trans_df[:, :Upgrade_Lvl])
    nonzero_storage_indices = findall(x -> x > 0, storage_df[:, :Storage_Energy])
    
    @constraint(model,
        old_gamma[a in nonzero_trans_indices],
        gamma[a] >= trans_df[a, :Upgrade_Lvl]
    )
    @constraint(model,
        old_s_energy[i in nonzero_storage_indices],
        s_energy_int[i] >= storage_df[i, :Storage_Energy]
    )
end

# Optimize method
import JuMP: optimize!

function optimize!(master::BendersMasterProblem)
    """
    Optimize the Benders master problem.
    """
    optimize!(master.jump_model)
end

# Get current investment values
function get_investments(master::BendersMasterProblem)
    """
    Get current investment decision values.
    Returns (gamma_val, s_energy_int_val)
    """
    gamma_val = value.(master.gamma)
    s_energy_int_val = value.(master.s_energy_int)
    return (gamma_val, s_energy_int_val)
end

# Get objective value
function get_objective_value(master::BendersMasterProblem)
    return objective_value(master.jump_model)
end

# Get theta values
function get_theta_values(master::BendersMasterProblem)
    return value.(master.theta)
end

# Add Benders cut
function add_benders_cut!(master::BendersMasterProblem, rep_index::Int, 
                         duals::Tuple, y_val::Tuple, phi_val::Float64)
    """
    Add a Benders optimality cut for a specific representative period.
    """
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_int_val = y_val
    
    @constraint(master.jump_model, 
        master.theta[rep_index] >= phi_val + 
        sum(dual_gamma[a] * (master.gamma[a] - gamma_val[a]) for a=1:master.E) +
        sum(dual_energy[i] * (master.s_energy_int[i] - s_energy_int_val[i]) for i=1:master.N)
    )
end

# Add feasibility cut
function add_feasibility_cut!(master::BendersMasterProblem, 
                             duals::Tuple, y_val::Tuple, ue_val::Float64)
    """
    Add a Benders feasibility cut.
    """
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_int_val = y_val
    
    @constraint(master.jump_model, 
        master.ue_sum >= ue_val + 
        sum(dual_gamma[a] * (master.gamma[a] - gamma_val[a]) for a=1:master.E) +
        sum(dual_energy[i] * (master.s_energy_int[i] - s_energy_int_val[i]) for i=1:master.N)
    )
end

# Trust region management
function add_trust_region!(master::BendersMasterProblem)
    """
    Add trust region constraints for stabilization.
    """
    if master.stabilization != "trust_region"
        return
    end
    
    if master.iter == 0
        # Initialize trust region at over-invested point
        gamma_core, s_energy_int_core = zeros(master.E), zeros(master.N)
        if master.warmstart
            gamma_core, s_energy_int_core = master.over_invested_point
        end
        y_core = [gamma_core, s_energy_int_core]
        push!(master.y_trust, y_core)
        master.jump_model.ext[:l1_radius] = [1]
        
        y_trust = master.y_trust[end]
        gamma_trust, s_energy_int_trust = y_trust
        
        abs_diff = master.jump_model[:abs_diff]
        trans_abs_diff = master.jump_model[:trans_abs_diff]
        
        @constraint(master.jump_model,
            trans_abs_diff_ub[a in 1:master.E],
            trans_abs_diff[a] >= master.gamma[a] - gamma_trust[a])
        @constraint(master.jump_model,
            trans_abs_diff_lb[a in 1:master.E],
            trans_abs_diff[a] >= gamma_trust[a] - master.gamma[a])
        
        @constraint(master.jump_model,
            abs_diff_ub[i in 1:master.N],
            abs_diff[i] >= master.s_energy_int[i] - s_energy_int_trust[i])
        @constraint(master.jump_model,
            abs_diff_lb[i in 1:master.N],
            abs_diff[i] >= s_energy_int_trust[i] - master.s_energy_int[i])
        
        @constraint(master.jump_model,
            trans_abs_diff_total,
            sum(trans_abs_diff[a] for a in 1:master.E) == 0)
        @constraint(master.jump_model,
            abs_diff_total,
            sum(abs_diff[i] for i in 1:master.N) == 0)
        return
    else
        remove_trust_region!(master)
    end
    
    y_trust = master.y_trust[end]
    gamma_trust, s_energy_int_trust = y_trust
    
    abs_diff = master.jump_model[:abs_diff]
    trans_abs_diff = master.jump_model[:trans_abs_diff]
    
    @constraint(master.jump_model,
        trans_abs_diff_ub[a in 1:master.E],
        trans_abs_diff[a] >= master.gamma[a] - gamma_trust[a])
    @constraint(master.jump_model,
        trans_abs_diff_lb[a in 1:master.E],
        trans_abs_diff[a] >= gamma_trust[a] - master.gamma[a])
    
    @constraint(master.jump_model,
        abs_diff_ub[i in 1:master.N],
        abs_diff[i] >= master.s_energy_int[i] - s_energy_int_trust[i])
    @constraint(master.jump_model,
        abs_diff_lb[i in 1:master.N],
        abs_diff[i] >= s_energy_int_trust[i] - master.s_energy_int[i])
    
    @constraint(master.jump_model,
        trans_abs_diff_total,
        sum(trans_abs_diff[a] for a in 1:master.E) <= master.jump_model.ext[:l1_radius][end])
    @constraint(master.jump_model,
        abs_diff_total,
        sum(abs_diff[i] for i in 1:master.N) == 2)
end

function remove_trust_region!(master::BendersMasterProblem)
    """
    Remove trust region constraints.
    """
    if master.stabilization == "trust_region"
        constraint_names = [
            :trans_abs_diff_ub, :trans_abs_diff_lb, :abs_diff_ub, :abs_diff_lb, 
            :abs_diff_total, :trans_abs_diff_total
        ]
    elseif master.stabilization == "boxstep"
        constraint_names = [
            :trans_box_ub, :trans_box_lb, :stor_box_lb, :stor_box_ub
        ]
    else
        return
    end
    
    obj_dict = object_dictionary(master.jump_model)
    
    for name in constraint_names
        if haskey(obj_dict, name)
            try
                constraint_obj = obj_dict[name]
                JuMP.delete(master.jump_model, constraint_obj)
                JuMP.unregister(master.jump_model, name)
                println("Removed constraint: $name")
            catch e
                println("Warning: Could not remove $name: $e")
                try
                    JuMP.unregister(master.jump_model, name)
                catch
                end
            end
        end
    end
end

# Level set management
function add_level_set!(master::BendersMasterProblem, current_obj::Float64)
    """
    Add level set constraint for convergence acceleration.
    """
    if !master.level_set
        return
    end
    remove_level_set!(master)
    
    obj_expr = master.jump_model.ext[:obj_expr]
    @constraint(master.jump_model, level_set, obj_expr <= current_obj)
end

function remove_level_set!(master::BendersMasterProblem)
    """
    Remove level set constraint.
    """
    obj_dict = object_dictionary(master.jump_model)
    name = :level_set
    if haskey(obj_dict, name)
        try
            constraint_obj = obj_dict[name]
            JuMP.delete(master.jump_model, constraint_obj)
            JuMP.unregister(master.jump_model, name)
            println("Removed constraint: $name")
        catch e
            println("Warning: Could not remove $name: $e")
            try
                JuMP.unregister(master.jump_model, name)
            catch
            end
        end
    end
end

# Iteration management
function increment_iteration!(master::BendersMasterProblem)
    """
    Increment the iteration counter.
    """
    master.iter += 1
end

function update_tracking!(master::BendersMasterProblem, y_val::Tuple, 
                         total_ue::Float64, total_obj::Float64)
    """
    Update tracking vectors for convergence monitoring.
    """
    master.last_y_val = [y_val[1], y_val[2]]
    push!(master.total_ue, total_ue)
    push!(master.total_obj, total_obj)
end