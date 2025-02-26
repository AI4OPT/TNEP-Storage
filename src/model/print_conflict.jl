"""
    print_conflict!(model)
Compute and print a conflict for an infeasible `model`.
"""
function print_conflict!(model)
    JuMP.compute_conflict!(model)
    ctypes = list_of_constraint_types(model)
    for (F, S) in ctypes
        cons = all_constraints(model, F, S)
        for i in eachindex(cons)
            isassigned(cons, i) || continue
            con = cons[i]
            cst = MOI.get(model, MOI.ConstraintConflictStatus(), con)
            cst == MOI.IN_CONFLICT && @info name(con) con
        end
    end
    return nothing
end


# Run these to identify between infeasible and unbounded
# set_optimizer_attribute(model, "InfUnbdInfo", 1)
# set_optimizer_attribute(model, "Presolve", 0)
# optimize!(model)

#=
If model is infeasible or unbounded
* Set InfUnbdInfo to 1, and re-optimize
* If that fails to give a specific status, disable presolve and re-optimize
* Once the model is known to be infeasible, run print_conflict!
=#

# function to print binding constraints
function print_binding_constraints_mip(model; filter_substring=nothing, output_file=nothing)
    tol = 1e-6
    binding_count = 0
    filtered_count = 0
    
    # Open file if specified
    io = isnothing(output_file) ? stdout : open(output_file, "w")
    
    println(io, "Binding constraints:")
    
    for (F, S) in list_of_constraint_types(model)
        for con in all_constraints(model, F, S)
            if isnothing(name(con)) || name(con) == ""
                continue
            end
            
            if !isnothing(filter_substring) && !occursin(filter_substring, name(con))
                continue
            end
            
            con_obj = constraint_object(con)
            
            try
                lhs_val = value(con)
                
                if S <: MOI.LessThan
                    rhs_val = MOI.constant(con_obj.set)
                    slack = rhs_val - lhs_val
                    sense = "<="
                elseif S <: MOI.GreaterThan
                    rhs_val = MOI.constant(con_obj.set)
                    slack = lhs_val - rhs_val
                    sense = ">="
                elseif S <: MOI.EqualTo
                    rhs_val = MOI.constant(con_obj.set)
                    slack = abs(lhs_val - rhs_val)
                    sense = "=="
                else
                    continue
                end
                
                if abs(slack) < tol
                    binding_count += 1
                    
                    if isnothing(filter_substring) || occursin(filter_substring, name(con))
                        println(io, "  $(name(con)): $lhs_val $sense $rhs_val (slack: $slack)")
                        filtered_count += 1
                    end
                end
            catch e
                continue
            end
        end
    end
    
    if isnothing(filter_substring)
        println(io, "Found $binding_count binding constraints")
    else
        println(io, "Found $filtered_count binding constraints containing '$filter_substring' (out of $binding_count total)")
    end
    
    # Close file if we opened one
    if !isnothing(output_file)
        close(io)
    end
    
    return nothing
end