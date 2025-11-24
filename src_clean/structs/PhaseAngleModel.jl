using JuMP
using CSV, DataFrames
include("OptimizationModel.jl")

# Phase Angle Model with direct flow constraints
mutable struct PhaseAngleModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String
    
    # Phase angle-specific fields
    solve_time::Float64
    
    # Pre-computed mappings for efficiency
    gen_bus_map::Dict{Int, Int}
    fbus::Vector{Int}
    tbus::Vector{Int}
    br_x::Vector{Float64}
    
    # Dimension sizes for convenience
    R::Int  # num_representatives
    N::Int  # num buses
    E::Int  # num branches
    T::Int  # num hours
    G::Int  # num generators
    K::Int  # num_cap_upgrades_max
    
    # Inner constructor
    function PhaseAngleModel(data::Dict{String, Any}, optimizer, simdir::String)
        
        # Create the base JuMP model
        jump_model = create_base_model(data, optimizer)
        
        # Get dimension sizes
        R = data["param"]["num_representatives"]
        N = length(data["bus"])
        E = length(data["branch"])
        T = data["param"]["num_hours"]
        G = length(data["gen"])
        K = data["param"]["num_cap_upgrades_max"]
        
        # Pre-compute gen_bus_map
        gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] 
                          for g in keys(data["gen"]))
        
        # Pre-compute branch properties
        fbus = [data["branch"]["$a"]["f_bus"] for a in 1:E]
        tbus = [data["branch"]["$a"]["t_bus"] for a in 1:E]
        br_x = [data["branch"]["$a"]["br_x"] for a in 1:E]
        
        # Add phase angle specific variables
        @variable(jump_model, pf[r in 1:R, a in 1:E, t in 1:T])
        @variable(jump_model, va[r in 1:R, i in 1:N, t in 1:T])
        
        # Add phase angle specific constraints
        add_phase_angle_constraints!(jump_model, data, R, N, E, T, fbus, tbus, br_x)
        
        # Create and return the model
        new(jump_model, data, simdir, 0.0, gen_bus_map, fbus, tbus, br_x, 
            R, N, E, T, G, K)
    end
end

# Add phase angle specific constraints
function add_phase_angle_constraints!(model::JuMP.Model, data::Dict{String, Any}, 
                                     R::Int, N::Int, E::Int, T::Int,
                                     fbus::Vector{Int}, tbus::Vector{Int}, 
                                     br_x::Vector{Float64})
    
    # Get rate_a nonzero branches
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    rate_a_zero = Set(parse(Int, x) for x in rate_a_zero)
    rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero)
    
    # Flow constraints
    @constraint(model, 
        flow_lb[r in 1:R, a in rate_a_nonzero, t in 1:T],
        model[:pf][r,a,t] >= -1 * data["branch"]["$a"]["rate_a"] - 
                             model[:gamma][a] * get_capacity_increment(data, a)
    )
    @constraint(model, 
        flow_ub[r in 1:R, a in rate_a_nonzero, t in 1:T],
        model[:pf][r,a,t] <= data["branch"]["$a"]["rate_a"] + 
                             model[:gamma][a] * get_capacity_increment(data, a)
    )
    
    # If rate a is zero (unlimited), don't allow upgrades
    @constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        model[:gamma][a] == 0
    )
    
    # Ohm's law constraint
    @constraint(model,
        ohms_law[r in 1:R, a in 1:E, t in 1:T],
        model[:pf][r,a,t] == (model[:va][r,tbus[a],t] - model[:va][r,fbus[a],t]) / br_x[a]
    )
    
    # Slack bus voltage angle
    @constraint(model,
        slack_bus_voltage[r in 1:R, t in 1:T],
        model[:va][r, 1, t] == 0
    )
    
    # Nodal power balance constraints
    @constraint(model, 
        power_balance[r in 1:R, i in 1:N, t in 1:T],
        sum(model[:pg][r,g,t] for g in [num for gen_type in values(data["bus"]["$i"]["gen"]) 
                                        if isa(gen_type, Array) for num in gen_type]) +
        sum(model[:pf][r,a,t] for a in [a for a in data["arcs_from"]["$i"] 
                                        if data["branch"]["$a"]["t_bus"] == i]) + # inflow
        model[:dis][r,i,t] +
        model[:ue][r,i,t] ==
        sum(model[:pf][r,a,t] for a in [a for a in data["arcs_from"]["$i"] 
                                        if data["branch"]["$a"]["f_bus"] == i]) + # outflow
        data["bus"]["$i"]["load"]["$r"][t] +
        model[:oe][r,i,t] +
        model[:ch][r,i,t]
    )
end

# Main solve function
function solve!(model::PhaseAngleModel; configure_optimizer::Bool=true)
    # Configure optimizer if requested
    if configure_optimizer
        set_optimizer_attribute(model.jump_model, "LogFile", 
                               joinpath(model.simdir, "gurobi_logfile.log"))
        set_optimizer_attribute(model.jump_model, "MIPGap", 
                               model.data["param"]["mip_gap"])
    end
    
    # Solve model
    t0 = time()
    JuMP.optimize!(model.jump_model)
    model.solve_time = time() - t0
    
    # Save results
    save_solve_time(model.simdir, model.solve_time)
    save_power_injections(model.simdir, model.jump_model, model.data)
    
    # Check termination status
    st = termination_status(model.jump_model)
    if st ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        @warn "Model did not solve to optimality. Status: $st"
    end
    
    return model.jump_model
end

# Optional: Function to plug in power injections from previous solve
function plug_in_power_injections!(model::PhaseAngleModel)
    csv_file = joinpath(model.simdir, "power_injections.csv")
    if !isfile(csv_file)
        return
    end

    df = CSV.read(csv_file, DataFrame)
    
    # Override with values from df
    for row in eachrow(df)
        if row.variable == "pg"
            fix(model.jump_model[:pg][row.rep, row.gen, row.time], row.value; force=true)
        elseif row.variable == "ue"
            fix(model.jump_model[:ue][row.rep, row.bus, row.time], row.value; force=true)
        elseif row.variable == "ch"
            if row.value >= 0
                fix(model.jump_model[:ch][row.rep, row.bus, row.time], row.value; force=true)
            else
                fix(model.jump_model[:dis][row.rep, row.bus, row.time], abs(row.value); force=true)
            end
        end
    end
end

# Convenience constructor that matches old function signature
function create_model_phase_angle(simdir, data::Dict{String, Any}, optimizer)
    model = PhaseAngleModel(data, optimizer, simdir)
    solve!(model)
    return model.jump_model
end