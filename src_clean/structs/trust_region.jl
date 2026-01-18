# Trust region management with incremental changes for multistage optimization
function add_trust_region!(master)
    """
    Add trust region constraints for stabilization.
    """
    if master.stabilization != "trust_region"
        return
    end
    
    years = sort(master.years)
    
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
        
        # Initialize L1 radius to 0 for first iteration (forces over-invested point)
        master.jump_model.ext[:l1_radius] = [0.0]
        
        # Add constraints
        tr_mode = get(master.data["param"], "tr_mode", "incremental")
        if tr_mode == "absolute"
            add_trust_region_constraints_abs!(master, y_core, years)
        else
            add_trust_region_constraints_incr!(master, y_core, years)
        end
        
        return
    else
        remove_trust_region!(master)
        # Add constraints with updated trust center
        y_trust = master.y_trust[end]

        tr_mode = get(master.data["param"], "tr_mode", "incremental")
        if tr_mode == "absolute"
            add_trust_region_constraints_abs!(master, y_trust, years)
        elseif tr_mode == "incremental"
            add_trust_region_constraints_incr!(master, y_trust, years)
        end
    end
end

function add_trust_region_constraints_abs!(master, y_trust, years)
    """
    Helper function to add trust region constraints.
    Handles absolute constraints for ALL years.
    
    Transmission lines use radius from l1_radius[end]
    Storage uses constant radius of 2 (except first iteration where it's 0)
    """
    println("[DEBUG] Adding ABSOLUTE trust regions for each period")
    abs_diff = master.jump_model[:abs_diff]
    trans_abs_diff = master.jump_model[:trans_abs_diff]
    trans_radius = master.jump_model.ext[:l1_radius][end]
    storage_radius = master.iter == 1 ? 0.0 : 2.0
    
    # Apply absolute trust region to ALL years: |x_y - x̌_y| <= r
    @constraint(master.jump_model,
        trans_abs_diff_ub[a in 1:master.E, year in years],
        trans_abs_diff[a, year] >= master.gamma[a, year] - y_trust[year][1][a])
    @constraint(master.jump_model,
        trans_abs_diff_lb[a in 1:master.E, year in years],
        trans_abs_diff[a, year] >= y_trust[year][1][a] - master.gamma[a, year])
    @constraint(master.jump_model,
        abs_diff_ub[i in 1:master.N, year in years],
        abs_diff[i, year] >= master.s_energy[i, year] - y_trust[year][2][i])
    @constraint(master.jump_model,
        abs_diff_lb[i in 1:master.N, year in years],
        abs_diff[i, year] >= y_trust[year][2][i] - master.s_energy[i, year])
    @constraint(master.jump_model,
        trans_abs_diff_total[year in years],
        sum(trans_abs_diff[a, year] for a in 1:master.E) <= trans_radius)
    @constraint(master.jump_model,
        abs_diff_total[year in years],
        sum(abs_diff[i, year] for i in 1:master.N) <= storage_radius)
end

function add_trust_region_constraints_incr!(master, y_trust, years)
    """
    Helper function to add trust region constraints.
    Handles both absolute (first year) and incremental (subsequent years) constraints.
    
    Transmission lines use radius from l1_radius[end]
    Storage uses constant radius of 2 (except first iteration where it's 0)
    """
    println("[DEBUG] Adding trust regions for each period by INCREMENTAL change")
    abs_diff = master.jump_model[:abs_diff]
    trans_abs_diff = master.jump_model[:trans_abs_diff]
    trans_radius = master.jump_model.ext[:l1_radius][end]
    storage_radius = master.iter == 1 ? 0.0 : 2.0
    
    first_year = years[1]
    
    # First year: absolute trust region |x_1 - x̌_1| <= r
    @constraint(master.jump_model,
        trans_abs_diff_ub_y1[a in 1:master.E],
        trans_abs_diff[a, first_year] >= master.gamma[a, first_year] - y_trust[first_year][1][a])
    @constraint(master.jump_model,
        trans_abs_diff_lb_y1[a in 1:master.E],
        trans_abs_diff[a, first_year] >= y_trust[first_year][1][a] - master.gamma[a, first_year])
    @constraint(master.jump_model,
        abs_diff_ub_y1[i in 1:master.N],
        abs_diff[i, first_year] >= master.s_energy[i, first_year] - y_trust[first_year][2][i])
    @constraint(master.jump_model,
        abs_diff_lb_y1[i in 1:master.N],
        abs_diff[i, first_year] >= y_trust[first_year][2][i] - master.s_energy[i, first_year])
    @constraint(master.jump_model,
        trans_abs_diff_total_y1,
        sum(trans_abs_diff[a, first_year] for a in 1:master.E) <= trans_radius)
    @constraint(master.jump_model,
        abs_diff_total_y1,
        sum(abs_diff[i, first_year] for i in 1:master.N) <= storage_radius)
    
    # Subsequent years: incremental trust region |(x_y - x_{y-1}) - (x̌_y - x̌_{y-1})| <= r
    if length(years) > 1
        K = 2:length(years)  # indices for "current year" positions

        # Incremental transmission upgrades:
        # trans_abs_diff[a, y_curr] >=  (γ_y - γ_{y-1}) - (γ̌_y - γ̌_{y-1})
        @constraint(master.jump_model,
            trans_abs_diff_incr_ub[a in 1:master.E, k in K],
            trans_abs_diff[a, years[k]] >=
                (master.gamma[a, years[k]] - master.gamma[a, years[k-1]]) -
                (y_trust[years[k]][1][a] - y_trust[years[k-1]][1][a])
        )

        # trans_abs_diff[a, y_curr] >= -( (γ_y - γ_{y-1}) - (γ̌_y - γ̌_{y-1}) )
        @constraint(master.jump_model,
            trans_abs_diff_incr_lb[a in 1:master.E, k in K],
            trans_abs_diff[a, years[k]] >=
                (y_trust[years[k]][1][a] - y_trust[years[k-1]][1][a]) -
                (master.gamma[a, years[k]] - master.gamma[a, years[k-1]])
        )

        # Incremental storage energy:
        # abs_diff[i, y_curr] >= (s_y - s_{y-1}) - (š_y - š_{y-1})
        @constraint(master.jump_model,
            abs_diff_incr_ub[i in 1:master.N, k in K],
            abs_diff[i, years[k]] >=
                (master.s_energy[i, years[k]] - master.s_energy[i, years[k-1]]) -
                (y_trust[years[k]][2][i] - y_trust[years[k-1]][2][i])
        )

        # abs_diff[i, y_curr] >= -( (s_y - s_{y-1}) - (š_y - š_{y-1}) )
        @constraint(master.jump_model,
            abs_diff_incr_lb[i in 1:master.N, k in K],
            abs_diff[i, years[k]] >=
                (y_trust[years[k]][2][i] - y_trust[years[k-1]][2][i]) -
                (master.s_energy[i, years[k]] - master.s_energy[i, years[k-1]])
        )

        # L1 norm constraints for incremental changes (per current year)
        @constraint(master.jump_model,
            trans_abs_diff_total_incr[k in K],
            sum(trans_abs_diff[a, years[k]] for a in 1:master.E) <= trans_radius
        )

        @constraint(master.jump_model,
            abs_diff_total_incr[k in K],
            sum(abs_diff[i, years[k]] for i in 1:master.N) <= storage_radius
        )
    end
end

function remove_trust_region!(master)
    """
    Remove trust region constraints for all years.
    Handles both absolute and incremental constraint modes.
    """
    if master.stabilization != "trust_region"
        return
    end

    tr_mode = get(master.data["param"], "tr_mode", "incremental")
    
    # Define constraint names based on mode
    if tr_mode == "absolute"
        constraint_names = [
            # Absolute constraints (applied to all years)
            :trans_abs_diff_ub, :trans_abs_diff_lb, 
            :abs_diff_ub, :abs_diff_lb,
            :abs_diff_total, :trans_abs_diff_total
        ]
    else  # incremental mode
        constraint_names = [
            # First year absolute constraints
            :trans_abs_diff_ub_y1, :trans_abs_diff_lb_y1, 
            :abs_diff_ub_y1, :abs_diff_lb_y1,
            :abs_diff_total_y1, :trans_abs_diff_total_y1,
            # Incremental constraints for subsequent years
            :trans_abs_diff_incr_ub, :trans_abs_diff_incr_lb,
            :abs_diff_incr_ub, :abs_diff_incr_lb,
            :abs_diff_total_incr, :trans_abs_diff_total_incr
        ]
    end
    
    obj_dict = object_dictionary(master.jump_model)
    
    for name in constraint_names
        if haskey(obj_dict, name)
            try
                constraint_obj = obj_dict[name]
                # Handle both single constraints and containers
                if constraint_obj isa JuMP.ConstraintRef
                    JuMP.delete(master.jump_model, constraint_obj)
                else
                    # It's a container - delete all constraints in it
                    JuMP.delete.(master.jump_model, constraint_obj)
                end
                JuMP.unregister(master.jump_model, name)
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