using SparseArrays
using LinearAlgebra

REF_IDX = 1

function create_susceptance_admittance_matrix(data::Dict{String, Any})
    branches = data["branch"]
    num_branches = length(branches)

    # Read branch info
    f_bus = [branches[string(i)]["f_bus"] for i in 1:num_branches]
    t_bus = [branches[string(i)]["t_bus"] for i in 1:num_branches]
    x = [branches[string(i)]["br_x"] for i in 1:num_branches]
    inv_x = [1/reactance for reactance in x]
    # inv_x = min.(inv_x, 500)

    # create diag(B) matrix
    diag_B = spdiagm(0 => inv_x)

    # create admittance matrix A
    num_buses = length(data["bus"]) 

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    # Populate the entries of the admittance matrix
    for e in 1:num_branches
        push!(rows, e)          # Branch index (row)
        push!(cols, f_bus[e])   # From-bus index (col)
        push!(vals, 1.0)        # +1 for from-bus

        push!(rows, e)          # Branch index (row)
        push!(cols, t_bus[e])   # To-bus index (col)
        push!(vals, -1.0)       # -1 for to-bus
    end
    A = sparse(rows, cols, vals, num_branches, num_buses)

    return diag_B, A
end

function compute_ptdf(diag_B, A)
    # Step 1: Compute A_T * diag_B * A
    L = A' * diag_B * A  # N x N

    # Adjust row and col corr. to slack bus
    L_reduced = L[2:end, 2:end]

    # Step 2: Invert the matrix...
    F = ldlt(L_reduced)  
    
    # Step 3: Solve the system
    dim = size(A, 2) - 1
    L_inv = Matrix{Float64}(I, dim, dim)  # Identity matrix

    for i in 1:dim
        L_inv[:, i] = F \ L_inv[:, i]
    end

    # Add back the slack bus dimension
    new_size = size(L_inv, 1) + 1
    L_inv_expanded = zeros(eltype(L_inv), new_size, new_size)
    L_inv_expanded[2:end, 2:end] = L_inv

    # Step 4: Get PTDF
    PTDF = diag_B * A * L_inv_expanded  # E x N

    return PTDF
end

function do_all_ptdf(data)
    diag_B, A_r = create_susceptance_admittance_matrix(data)
    ptdf = compute_ptdf(diag_B, A_r)

    return ptdf
end

function do_all_incidence(data)
    branches = data["branch"]
    buses = data["bus"]

    num_branches = length(branches)
    num_buses = length(buses)

    f_bus = [branches[string(i)]["f_bus"] for i in 1:num_branches]
    t_bus = [branches[string(i)]["t_bus"] for i in 1:num_branches]

    # Prepare to build incidence matrix (rows = buses, cols = branches)
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for (l, (f, t)) in enumerate(zip(f_bus, t_bus))
        push!(rows, f)   # From bus: -1
        push!(cols, l)
        push!(vals, -1.0)

        push!(rows, t)   # To bus: +1
        push!(cols, l)
        push!(vals, +1.0)
    end

    A = sparse(rows, cols, vals, num_buses, num_branches)
    return A
end
