# Abstract base type for all models
abstract type OptimizationModel end

# Specific model types
struct PhaseAngleModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String
end

struct PTDFModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String
    
    #PTDF-specific fields
    ptdf_matrix::Matrix{Float64}
    tracked_constraints::Dict{Tuple{Int,Int,Int,Bool}, Bool}
    rate_a_nonzero::Set{Int}
    max_ptdf_iterations::Int
    max_ptdf_per_iteration::Int
    ptdf_tol::Float64
    solve_time::Float64
    
    function PTDFModel(jump_model, data, simdir; 
                                max_ptdf_iterations=256,
                                max_ptdf_per_iteration=32,
                                ptdf_tol=1e-6)

        ptdf_matrix = do_all_ptdf(data)
        
        if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
            ptdf_sparse = map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"]
            ptdf_matrix = ptdf_matrix .* ptdf_sparse
        end
        
        rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
        rate_a_set = Set(parse(Int, x) for x in rate_a_nonzero)
        
        new(jump_model, data, simdir, ptdf_matrix, 
            Dict{Tuple{Int,Int,Int,Bool}, Bool}(),
            rate_a_set, max_ptdf_iterations, max_ptdf_per_iteration, 
            ptdf_tol, 0.0)
    end
end