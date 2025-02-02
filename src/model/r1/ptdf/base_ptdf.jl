using SparseArrays
using LinearAlgebra

function create_susceptance_admittance_matrix(data::Dict{String, Any})
    branches = data["branch"]
    num_branches = length(branches)

    # Read branch info
    f_bus = [branches[string(i)]["f_bus"] for i in 1:num_branches]
    t_bus = [branches[string(i)]["t_bus"] for i in 1:num_branches]
    x = [branches[string(i)]["br_x"] for i in 1:num_branches]
    inv_x = [1/reactance for reactance in x]

    # create diag(B) matrix
    diag_B = spdiagm(0 => inv_x)

    # create reduced admittance matrix A_r
    num_buses = length(data["bus"]) 

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    # Populate the entries of the admittance matrix
    for e in 1:num_branches
        if f_bus[e] != 1  # Exclude the slack bus
            push!(rows, e)          # Branch index (row)
            push!(cols, f_bus[e]-1) # Adjust for N-1 (exclude slack bus)
            push!(vals, 1.0)        # +1 for from-bus
        end

        if t_bus[e] != 1  # Exclude the slack bus
            push!(rows, e)          # Branch index (row)
            push!(cols, t_bus[e]-1) # Adjust for N-1 (exclude slack bus)
            push!(vals, -1.0)       # -1 for to-bus
        end
    end
    A_r = sparse(rows, cols, vals, num_branches, num_buses-1)

    return diag_B, A_r
end

function compute_ptdf(diag_B, A_r)
    # Step 1: Compute A_r_T * diag_B * A_r
    L = A_r' * diag_B * A_r  # (N-1) x (N-1)
    # Step 2: Invert the matrix
    L_inv = inv(Matrix(L))  # Convert to dense for inversion
    # Step 3: Compute PTDF
    PTDF = diag_B * A_r * L_inv  # E x (N-1)

    return PTDF
end

function do_all_ptdf(data)
    diag_B, A_r = create_susceptance_admittance_matrix(data)
    ptdf = compute_ptdf(diag_B, A_r)

    return ptdf
end
