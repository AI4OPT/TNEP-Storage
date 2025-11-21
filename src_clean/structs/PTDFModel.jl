using JuMP
using CSV, DataFrames
include("OptimizationModel.jl")

# PTDF Model with iterative constraint addition
mutable struct PTDFModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String
    
    # PTDF-specific fields
    ptdf_matrix::Matrix{Float64}
    tracked_constraints::Dict{Tuple{Int,Int,Int,Bool}, Bool}
    rate_a_nonzero::Set{Int}
    max_ptdf_iterations::Int
    max_ptdf_per_iteration::Int
    ptdf_tol::Float64
    solve_time::Float64
    
    # Pre-computed mappings for efficiency
    gen_bus_map::Dict{Int, Int}
    
    # Dimension sizes for convenience
    R::Int  # num_representatives
    N::Int  # num buses
    E::Int  # num branches
    T::Int  # num hours
    G::Int  # num generators
    
    # Inner constructor
    function PTDFModel(data::Dict{String, Any}, optimizer, simdir::String; 
                      max_ptdf_iterations::Int=256,
                      max_ptdf_per_iteration::Int=32,
                      ptdf_tol::Float64=1e-6)
        
        # Create the base JuMP model
        jump_model = create_base_model(data, optimizer)
        
        # Compute PTDF matrix
        ptdf_matrix = do_all_ptdf(data)
        
        # Apply cutoff if specified
        if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
            ptdf_sparse = map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"]
            ptdf_matrix = ptdf_matrix .* ptdf_sparse
        end
        
        # Get rate_a nonzero branches
        rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
        rate_a_set = Set(parse(Int, x) for x in rate_a_nonzero)
        
        # Pre-compute gen_bus_map
        gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] 
                          for g in keys(data["gen"]))
        
        # Get dimension sizes
        R = data["param"]["num_representatives"]
        N = length(data["bus"])
        E = length(data["branch"])
        T = data["param"]["num_hours"]
        G = length(data["gen"])
        
        # Add container for flow constraints to JuMP model
        jump_model[:ptdf_flow] = Dict{String, ConstraintRef}()
        
        # Create and return the model
        new(jump_model, data, simdir, ptdf_matrix, 
            Dict{Tuple{Int,Int,Int,Bool}, Bool}(),
            rate_a_set, max_ptdf_iterations, max_ptdf_per_iteration, 
            ptdf_tol, 0.0, gen_bus_map, R, N, E, T, G)
    end
end

# Warm-start from saved constraints
function warm_start!(model::PTDFModel)
    constraint_file = joinpath(model.simdir, "tracked_constraints.csv")
    if !isfile(constraint_file)
        return
    end
    
    tracked_df = CSV.read(constraint_file, DataFrame)
    println("Beginning optimized warm-start process")
    
    # Group the constraints by arc once
    arc_groups = groupby(tracked_df, :arc)
    total_constraints = 0
    
    for group in arc_groups
        a = first(group).arc
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        cap_increment = get_capacity_increment(model.data, a)
        ptdf_row = model.ptdf_matrix[a, :]
        
        # Extract all the (r,t) pairs for this arc
        rep_time_pairs = [(cons.rep, cons.time, cons.ub) for cons in eachrow(group)]
        
        # Pre-compute the line limits for all constraints in this group
        line_limit = line_limit_base + model.jump_model[:gamma][a] * cap_increment
        
        # Create constraints for each (r,t) pair
        for (r, t, ub) in rep_time_pairs
            # Build the flow expression
            flow_expr = sum(ptdf_row[i] * (
                sum(model.jump_model[:pg][r, g, t] for g in 1:model.G 
                    if model.gen_bus_map[g] == i; init=0.0) -
                model.data["bus"]["$i"]["load"]["$r"][t] +
                model.jump_model[:ue][r, i, t] -
                model.jump_model[:ch][r, i, t] +
                model.jump_model[:dis][r, i, t]
            ) for i in 1:model.N)
            
            # Add constraint
            if ub
                model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = 
                    @constraint(model.jump_model, flow_expr <= line_limit)
            else
                model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = 
                    @constraint(model.jump_model, flow_expr >= -line_limit)
            end
            
            # Track the constraint
            model.tracked_constraints[(a, r, t, ub)] = true
            total_constraints += 1
        end
    end
    
    println("Warm-started with $total_constraints constraints")
end

# Find violations in current solution
function find_violations(model::PTDFModel, flows, gamma_values)
    violations = []
    
    for a in model.rate_a_nonzero
        # Pre-compute for this branch
        cap_increment = get_capacity_increment(model.data, a)
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        line_limit_fixed = line_limit_base + gamma_values[a] * cap_increment
        
        for r in 1:model.R, t in 1:model.T
            flow_val = flows[a, r, t]
            violation_amount = abs(flow_val) - line_limit_fixed
            ub = (flow_val >= 0)

            # Skip if already tracked
            if get(model.tracked_constraints, (a, r, t, ub), false)
                continue
            end
    
            if violation_amount > model.ptdf_tol
                push!(violations, (a, r, t, ub, violation_amount))
            end
        end
    end
    
    return violations
end

# Add PTDF constraints for violations
function add_constraints!(model::PTDFModel, sorted_violations)
    n_added = 0
    max_to_add = model.max_ptdf_per_iteration
    
    for (a, r, t, ub, _) in Iterators.take(sorted_violations, max_to_add)
        ptdf_row = model.ptdf_matrix[a, :]
        cap_increment = get_capacity_increment(model.data, a)
        line_limit_base = model.data["branch"]["$a"]["rate_a"]
        line_limit_expr = line_limit_base + model.jump_model[:gamma][a] * cap_increment
    
        # Define the flow expression
        flow_expr = sum(ptdf_row[i] * (
            sum(model.jump_model[:pg][r, g, t] for g in 1:model.G 
                if model.gen_bus_map[g] == i; init=0.0)
            - model.data["bus"]["$i"]["load"]["$r"][t]
            + model.jump_model[:ue][r, i, t] 
            - model.jump_model[:ch][r, i, t] 
            + model.jump_model[:dis][r, i, t]
        ) for i in 1:model.N)
    
        # Add constraint
        if ub
            model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = 
                @constraint(model.jump_model, flow_expr <= line_limit_expr)
        else
            model.jump_model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = 
                @constraint(model.jump_model, flow_expr >= -line_limit_expr)
        end
    
        model.tracked_constraints[(a, r, t, ub)] = true
        n_added += 1
    end
    
    return n_added
end

# Main solve function with iterative constraint addition
function solve!(model::PTDFModel; configure_optimizer::Bool=false)
    # Configure optimizer if requested
    if configure_optimizer
        set_optimizer_attribute(model.jump_model, "LogFile", 
                               joinpath(model.simdir, "gurobi_logfile.log"))
        set_optimizer_attribute(model.jump_model, "MIPGap", 
                               model.data["param"]["mip_gap"])
    end
    
    # Warm-start from previous constraints if available
    warm_start!(model)
    
    # Begin lazy PTDF loop
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < model.max_ptdf_iterations
        # Solve model
        JuMP.optimize!(model.jump_model)
        
        # Exit if not solved optimally
        st = termination_status(model.jump_model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(model.jump_model[:pg])
        ue_values = value.(model.jump_model[:ue])
        ch_values = value.(model.jump_model[:ch])
        dis_values = value.(model.jump_model[:dis])
        gamma_values = value.(model.jump_model[:gamma])
        
        # Compute all flows at once
        flows = zeros(length(model.rate_a_nonzero), model.R, model.T)
        compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, 
                      model.data, model.ptdf_matrix)

        # Find violations
        violations = find_violations(model, flows, gamma_values)
        sorted_violations = sort(violations, by = x -> -x[5])
        n_violated = length(sorted_violations)
        
        # Add constraints for top violations
        n_added = add_constraints!(model, sorted_violations)

        if n_added > 0
            save_tracked_constraints(model.simdir, model, n_violated)
        end
        
        solved = (n_violated == 0)
        niter += 1
        println("Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    model.solve_time = time() - t0
    save_solve_time(model.simdir, model.solve_time)
    save_power_injections(model.simdir, model.jump_model, model.data)
    
    return model.jump_model
end

# Helper function to compute flows
function compute_flows!(flows, pg_values, ue_values, ch_values, dis_values, data, PTDF)
    n_buses = length(data["bus"])
    net_injections = zeros(n_buses)
    
    for r in 1:size(pg_values, 1), t in 1:size(pg_values, 3)
        fill!(net_injections, 0.0)

        for i in 1:n_buses 
            net_injections[i] = sum(pg_values[r, g, t] 
                                   for g in 1:length(data["gen"]) 
                                   if data["gen"]["$g"]["gen_bus"] == i; 
                                   init=0.0)
            net_injections[i] -= data["bus"]["$i"]["load"]["$r"][t]
            net_injections[i] += ue_values[r, i, t]
            net_injections[i] -= ch_values[r, i, t]
            net_injections[i] += dis_values[r, i, t]
        end
        
        flows[:, r, t] = PTDF * net_injections
    end
end

# Convenience constructor that matches old function signature
function create_model_r1_ptdf_iterative_simplified_sorted(
    simdir, data::Dict{String, Any}, optimizer; 
    max_ptdf_iterations::Int=256,
    max_ptdf_per_iteration::Int=32,
    ptdf_tol::Float64=1e-6)
    
    model = PTDFModel(data, optimizer, simdir;
                     max_ptdf_iterations=max_ptdf_iterations,
                     max_ptdf_per_iteration=max_ptdf_per_iteration,
                     ptdf_tol=ptdf_tol)
    
    solve!(model)
    
    return model.jump_model
end