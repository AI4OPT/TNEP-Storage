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
set_optimizer_attribute(model, "InfUnbdInfo", 1)
set_optimizer_attribute(model, "Presolve", 0)
optimize!(model)

#=
If model is infeasible or unbounded
* Set InfUnbdInfo to 1, and re-optimize
* If that fails to give a specific status, disable presolve and re-optimize
* Once the model is known to be infeasible, run print_conflict!
=#