using JuMP
using JuMP.Containers
using Gurobi
using Serialization
using CSV, DataFrames

# Multistage Benders Master Problem struct
mutable struct MultistageBendersMasterProblem
    jump_model::JuMP.Model
    data::Dict{String, Any}
    superdir::String
    
    # Decision variables (stored for convenience)
    gamma::DenseAxisArray{VariableRef, 2}  # gamma[arc, year]
    s_energy::DenseAxisArray{VariableRef, 2}  # s_energy[bus, year]
    theta::DenseAxisArray{VariableRef, 1} # theta[string(date)]
    ue_sum::VariableRef
    
    # Metadata and tracking 
    date_weights::Dict{String, Float64}
    y_trust::Vector{Dict{Int, Vector{Vector{Float64}}}} # e.g. [Dict(2030 => [gamma_2030, s_energy_2030]), ...]
    gap::Vector{Float64}
    iter::Int
    total_ue::Vector{Float64}
    total_obj::Vector{Float64}
    upper_bound::Float64
    lower_bound::Float64
    last_y_val::Dict{Int, Vector{Vector{Float64}}} # e.g. Dict(2030 => [gamma_2030, s_energy_2030])
    over_invested_point::Dict{Int, Vector{Vector{Float64}}} # e.g. Dict(2030 => [gamma_2030, s_energy_2030])
    
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
    years::Vector{Int}  # list of years

    # Discount factors by year e.g. disc_factors[2035] = 0.78
    disc_factors::Dict{Int, Float64}()
    
    # Inner constructor
    function MultistageBendersMasterProblem(superdir::String, data::Dict{String, Any}, 
                                 date_weights::Dict{String, Float64})
        
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

        # Compute discount factor
        discount_rate = get(data["param"], "discount_rate", 0.0)
        disc_factor = 1 / (1 + discount_rate)
        
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

        years = get(data["param"], "years", [data["param"]["decarbonization_year"]])        
        # Investment level of capacity upgrade
        @variable(jump_model, 0 <= gamma[a=1:E, y in years] <= K, Int)
        
        # Energy rating of storage
        @variable(jump_model, s_energy[i=1:N, y in years] >= 0, Int)
        
        # Subproblem objectives
        dates = [string(year) * x[5:end] for year in years for x in data["param"]["dates"]]
        @variable(jump_model, theta[dates] >= 0)
 
        # Subproblem load shed
        @variable(jump_model, ue_sum >= 0)
        
        # Trust region variables if needed
        if stabilization == "trust_region"
            @variable(jump_model, abs_diff[1:N, y in years] >= 0)
            @variable(jump_model, trans_abs_diff[1:E, y in years] >= 0)
        end
        
        #
        #   II. Constraints
        #

        # Monotonicity constraints: investments can only increase over time
        if Y > 1
            years = sort(years)
            for year_idx in 2:length(years)
                y_curr = years[year_idx]
                y_prev = years[year_idx - 1]
                
                @constraint(jump_model,
                    monotone_gamma[a=1:E],
                    gamma[a, y_curr] >= gamma[a, y_prev]
                )
                
                @constraint(jump_model,
                    monotone_storage[i=1:N],
                    s_energy[i, y_curr] >= s_energy[i, y_prev]
                )
            end
        end
        
        if warmstart
            @constraint(jump_model, no_load_shed, ue_sum <= 0)
        end
        
        # Check if previous investments exist
        """
        if haskey(data["param"], "previous_investment_dir")
            prev_dir = data["param"]["previous_investment_dir"]
            add_prev_upgrades_internal!(superdir, jump_model, data, prev_dir, gamma, s_energy)
        end
        """
        
        # If rate a is zero (unlimited), then don't allow upgrades (for any year)
        @constraint(jump_model, 
            rate_a_zero_line_upgrade[a in rate_a_zero, y in years],
            gamma[a, y] == 0
        )
        
        # Energy rating maximum (for every year)
        @constraint(jump_model, 
            installed_energy_ub[i in 1:N, y in years],
            s_energy[i, y] * data["param"]["storage_energy_size"] <= data["param"]["max_energy_rating"]
        )
        
        # Compute over-invested point
        over_invested_point = compute_superset_core_point(superdir, is_multistage=true)
        warmstart_dir = joinpath(superdir, "warmstart")
        mkpath(warmstart_dir)
        for year in years
            export_investments_csv(data, over_invested_points[year][1], over_invested_point[year][2], output_dir=joinpath(warmstart_dir, string(year)))
        end

        # Storage candidates 
        if get(data["param"], "storage_cand", false)
            storage_upgrades = over_invested_points[years[end]][2]
            storage_non_indices = findall(x -> x == 0, storage_upgrades)
            
            @constraint(jump_model, 
                energy_cand[i in storage_non_indices, y in years],
                s_energy[i, y] == 0
            )
        end
        
        # Line candidates
        if get(data["param"], "line_cand", false)
            line_upgrades = over_invested_points[years[end]][1]
            line_non_indices = findall(x -> x == 0, line_upgrades)
            
            @constraint(jump_model, 
                line_cand[a in line_non_indices, y in years],
                gamma[a, y] == 0
            )
        end
        
        #
        #   III. Objective
        #
        years = sort(years)
        first_year = first(years)

        disc_factors = Dict{Int, Float64}()
        for year in years
            disc_factors[year] = disc_factor^(year - first_year)
        end

        # Investment costs over time with discounting
        # First stage: c^T x_1
        first_stage_cost = (
            sum(s_energy[i, first_year] * data["param"]["storage_energy_size"] * 
                data["param"]["bess_energy_cost"] for i in 1:N) + 
            sum(data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * 
                data["branch"]["$a"]["distance"] * gamma[a, first_year] for a in 1:E)
        )

        # Subsequent stages: sum over y=2..Y of d^(y-1) * c^T (x_y - x_{y-1})
        incremental_costs = AffExpr(0.0)
        if Y > 1
            for year_idx in 2:length(years)
                y_curr = years[year_idx]
                y_prev = years[year_idx - 1]
                d = disc_factors[y_curr]
                
                add_to_expression!(incremental_costs,
                    d * sum(s_energy[i, y_curr] - s_energy[i, y_prev] for i in 1:N) * 
                    data["param"]["storage_energy_size"] * data["param"]["bess_energy_cost"]
                )
                
                add_to_expression!(incremental_costs,
                    d * sum((gamma[a, y_curr] - gamma[a, y_prev]) * 
                    data["param"]["cap_upgrade_cost"] * data["param"]["cap_upgrade_increment"] * 
                    data["branch"]["$a"]["distance"] for a in 1:E)
                )
            end
        end

        # Operational costs: sum over years of dates of (theta_{1,s} + d*theta_{2,s} + d^2*...)
        operational_costs = sum(
            disc_factors[y] * theta[date, y] * date_weights[date[6:end]] 
            for y in years, date in dates
        )

        obj_expr = (
            first_stage_cost + 
            incremental_costs +
            operational_costs + 
            ue_sum * data["param"]["operational_weight"] * data["param"]["under_served_penalty"]
        )
        
        @objective(jump_model, Min, obj_expr)
        
        # Store objective expression for later use
        jump_model.ext[:obj_expr] = obj_expr

        # create last_y_val
        last_y_val = Dict(year => [zeros(E), zeros(N)] for year in years)
        
        # Create and return the struct
        new(jump_model, data, superdir,
            gamma, s_energy, theta, ue_sum,
            date_weights,
            Vector{Dict{Int, Vector{Vector{Float64}}}}(),  # y_trust
            Vector{Float64}(),  # gap
            1,  # iter
            [Inf],  # total_ue
            [Inf],  # total_obj
            Inf, # upper_bound
            0.0, # lower_bound
            last_y_val,  # last_y_val
            over_invested_point,
            warmstart, stabilization, level_set_flag,
            R, N, E, T, G, K, years, disc_factors)
    end
end

function solve!(master::MultistageBendersMasterProblem)
    """
    Optimize the Benders master problem.
    """
    optimize!(master.jump_model)
end

# Helper function to get investments for a specific year
function get_investments(master::MultistageBendersMasterProblem, year::Int)
    """
    Get investment decision values for a specific year.
    Returns (gamma_val, s_energy_val)
    """
    gamma_val = [value(master.gamma[a, year]) for a in 1:master.E]
    s_energy_val = [value(master.s_energy[i, year]) for i in 1:master.N]
    return (gamma_val, s_energy_val)
end

# Get objective value
function get_objective_value(master::MultistageBendersMasterProblem)
    return objective_value(master.jump_model)
end

# Get theta values
function get_theta_values(master::MultistageBendersMasterProblem)
    return value.(master.theta)
end

# Add Benders cut
function add_benders_cut!(master::MultistageBendersMasterProblem, date::String, year::Int,
                         duals::Tuple, y_val::Tuple, phi_val::Float64)
    """
    Add a Benders optimality cut for a specific representative period and year.
    """
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_val = y_val
    
    @constraint(master.jump_model, 
        master.theta[date] >= phi_val + 
        sum(dual_gamma[a] * (master.gamma[a, year] - gamma_val[a]) for a=1:master.E) +
        sum(dual_energy[i] * (master.s_energy[i, year] - s_energy_val[i]) for i=1:master.N)
    )
end

# Add feasibility cut
function add_feasibility_cut!(master::MultistageBendersMasterProblem, year::Int,
                             duals::Tuple, y_val::Tuple, ue_val::Float64)
    """
    Add a Benders feasibility cut for a specific year.
    
    Note: For multistage, feasibility cuts may need to be reformulated
    depending on how infeasibility is handled across years.
    """
    dual_gamma, dual_energy = duals
    gamma_val, s_energy_val = y_val
    
    @constraint(master.jump_model, 
        master.ue_sum >= ue_val + 
        sum(dual_gamma[a] * (master.gamma[a, year] - gamma_val[a]) for a=1:master.E) +
        sum(dual_energy[i] * (master.s_energy[i, year] - s_energy_val[i]) for i=1:master.N)
    )
end

# Trust region management # TODO year-dependent l1 radii updates potentially needed
function add_trust_region!(master::MultistageBendersMasterProblem)
    """
    Add trust region constraints for stabilization.
    Year-specific trust regions |x_y - x̌_y| <= r_y for each year
    """
    if master.stabilization != "trust_region"
        return
    end
    
    years = master.years
    
    if master.iter == 1
        # Initialize trust region at over-invested point for each year
        y_core = Dict{Int, Vector{Vector{Float64}}}()
        for year in years
            if master.warmstart
                gamma_core, s_energy_core = master.over_invested_point[year]
            else
                gamma_core, s_energy_core = zeros(master.E), zeros(master.N)
            end
            y_core[year] = [gamma_core, s_energy_core]
        end
        push!(master.y_trust, y_core)
        master.jump_model.ext[:l1_radius] = [1]
        
        # Add constraints for each year
        y_trust = master.y_trust[end]
        abs_diff = master.jump_model[:abs_diff]
        trans_abs_diff = master.jump_model[:trans_abs_diff]
        
        for year in years
            gamma_trust, s_energy_trust = y_trust[year]
            
            @constraint(master.jump_model,
                trans_abs_diff_ub[a in 1:master.E, y=year],
                trans_abs_diff[a, y] >= master.gamma[a, y] - gamma_trust[a])
            @constraint(master.jump_model,
                trans_abs_diff_lb[a in 1:master.E, y=year],
                trans_abs_diff[a, y] >= gamma_trust[a] - master.gamma[a, y])
            
            @constraint(master.jump_model,
                abs_diff_ub[i in 1:master.N, y=year],
                abs_diff[i, y] >= master.s_energy[i, y] - s_energy_trust[i])
            @constraint(master.jump_model,
                abs_diff_lb[i in 1:master.N, y=year],
                abs_diff[i, y] >= s_energy_trust[i] - master.s_energy[i, y])
            
            # Initialize with zero radius for first iteration
            @constraint(master.jump_model,
                trans_abs_diff_total[y=year],
                sum(trans_abs_diff[a, y] for a in 1:master.E) == 0)
            @constraint(master.jump_model,
                abs_diff_total[y=year],
                sum(abs_diff[i, y] for i in 1:master.N) == 0)
        end
        return
    else
        remove_trust_region!(master)
    end
    
    # Add constraints for each year with proper radius
    y_trust = master.y_trust[end]
    abs_diff = master.jump_model[:abs_diff]
    trans_abs_diff = master.jump_model[:trans_abs_diff]
    
    for year in years
        gamma_trust, s_energy_trust = y_trust[year]
        
        @constraint(master.jump_model,
            trans_abs_diff_ub[a in 1:master.E, y=year],
            trans_abs_diff[a, y] >= master.gamma[a, y] - gamma_trust[a])
        @constraint(master.jump_model,
            trans_abs_diff_lb[a in 1:master.E, y=year],
            trans_abs_diff[a, y] >= gamma_trust[a] - master.gamma[a, y])
        
        @constraint(master.jump_model,
            abs_diff_ub[i in 1:master.N, y=year],
            abs_diff[i, y] >= master.s_energy[i, y] - s_energy_trust[i])
        @constraint(master.jump_model,
            abs_diff_lb[i in 1:master.N, y=year],
            abs_diff[i, y] >= s_energy_trust[i] - master.s_energy[i, y])
        
        @constraint(master.jump_model,
            trans_abs_diff_total[y=year],
            sum(trans_abs_diff[a, y] for a in 1:master.E) <= master.jump_model.ext[:l1_radius][end])
        @constraint(master.jump_model,
            abs_diff_total[y=year],
            sum(abs_diff[i, y] for i in 1:master.N) <= master.jump_model.ext[:l1_radius][end])
    end
end

# Level set management
function add_level_set!(master::MultistageBendersMasterProblem, current_obj::Float64)
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

function remove_level_set!(master::MultistageBendersMasterProblem)
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
function increment_iteration!(master::MultistageBendersMasterProblem)
    """
    Increment the iteration counter.
    """
    master.iter += 1
end

function update_tracking!(master::MultistageBendersMasterProblem, y_val::Dict{Int, Vector{Vector{Float64}}}, 
                         total_ue::Float64, total_obj::Float64)
    """
    Update tracking vectors for convergence monitoring.
    Tracks investments and objectives for each year/stage
    """
    master.last_y_val = y_val
    push!(master.total_ue, total_ue)
    push!(master.total_obj, total_obj)
end