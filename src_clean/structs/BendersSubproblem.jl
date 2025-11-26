include("PTDFModel.jl")

mutable struct BendersSubproblem
    ptdf_model::PTDFModel
    current_investments::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    
    function BendersSubproblem(data::Dict{String, Any}, optimizer, simdir::String;
                              max_ptdf_iterations::Int=256,
                              max_ptdf_per_iteration::Int=32,
                              ptdf_tol::Float64=1e-6)
        
        ptdf_model = PTDFModel(data, optimizer, simdir;
                              max_ptdf_iterations=max_ptdf_iterations,
                              max_ptdf_per_iteration=max_ptdf_per_iteration,
                              ptdf_tol=ptdf_tol)
        
        new(ptdf_model, nothing)
    end
end

function fix_investments!(sub::BendersSubproblem, gamma_val::Vector, s_energy_val::Vector)
    """
    Fix investment decisions for this Benders subproblem.
    """
    fix_investments!(sub.ptdf_model, gamma_val, s_energy_val)
    sub.current_investments = (gamma_val, s_energy_val)
end

function solve!(sub::BendersSubproblem, objective_type::Symbol; configure_optimizer::Bool=false)
    """
    Solve the Benders subproblem with specified objective.
    """
    if isnothing(sub.current_investments)
        error("Must fix investments before solving subproblem")
    end
    
    # Set objective
    set_objective!(sub.ptdf_model, objective_type)
    
    # Solve
    solve!(sub.ptdf_model; configure_optimizer=configure_optimizer)
    
    return sub.ptdf_model.jump_model
end

function extract_duals(sub::BendersSubproblem)
    """
    Extract duals for Benders cut generation.
    """
    return extract_investment_duals(sub.ptdf_model)
end

function get_objective_value(sub::BendersSubproblem)
    return objective_value(sub.ptdf_model.jump_model)
end

function get_total_unserved_energy(sub::BendersSubproblem)
    ue = value.(sub.ptdf_model.jump_model[:ue])
    return sum(ue)
end