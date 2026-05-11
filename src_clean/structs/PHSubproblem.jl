using JuMP
@isdefined(OptimizationModel) || include("OptimizationModel.jl")
@isdefined(PTDFModel)         || include("PTDFModel.jl")
include("../helpers/helpers.jl")

# ---------------------------------------------------------------------------
# PHBlockModel
#
# One temporal block (T_b hours) for a single representative period r.
# Mirrors PTDFModel's fields so PTDF lazy-constraint dispatch works unchanged.
#
# Boundary conventions:
#   block_id == 1         : normal soc_start, no soc_end (free right boundary)
#   1 < block_id < n_blocks : free soc_in_b (left boundary), no soc_end (free right)
#   block_id == n_blocks  : free soc_in_b (left boundary), normal soc_end
# ---------------------------------------------------------------------------
mutable struct PHBlockModel <: OptimizationModel
    jump_model::JuMP.Model
    data::Dict{String, Any}
    simdir::String

    ptdf_matrix::Matrix{Float64}
    tracked_constraints::Dict{Tuple{Int,Int,Int,Bool}, Bool}
    rate_a_nonzero::Set{Int}
    max_ptdf_iterations::Int
    max_ptdf_per_iteration::Int
    ptdf_tol::Float64
    solve_time::Float64
    constraint_stats::Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}
    gen_bus_map::Dict{Int, Int}

    # Dimensions — T here is T_b (block size), not the full horizon
    R::Int; N::Int; E::Int; T::Int; G::Int

    block_id::Int   # 1:n_blocks
    n_blocks::Int
    r::Int          # representative index (for load profile lookup)
    t_start::Int    # global hour of first period in this block
end

# ---------------------------------------------------------------------------
# PHSubproblem
#
# Sibling to BendersSubproblem. Given a fixed master first-stage solution,
# decomposes the operational subproblem into n_blocks temporal chunks per
# representative period and solves them via ADMM. The coupling variables are
# SOC values at each of the n_blocks-1 block boundaries.
# ---------------------------------------------------------------------------
mutable struct PHSubproblem
    block_models::Matrix{PHBlockModel}  # [R, n_blocks]
    n_blocks::Int

    current_investments::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}
    storage_nodes::Vector{Int}         # nodes where s_energy_val > 0 (set by fix_investments!)

    # ADMM state — all arrays [N, R, n_blocks-1], one slice per boundary
    lambda::Array{Float64, 3}    # Lagrange multipliers
    rho::Float64                 # augmented Lagrangian penalty
    soc_bar::Array{Float64, 3}   # consensus SOC at each boundary
    soc_out::Array{Float64, 3}   # soc[i, T_b] extracted from block bd after solve
    soc_in::Array{Float64, 3}    # soc_in_b[i] extracted from block bd+1 after solve

    max_ph_iterations::Int
    ph_tol::Float64

    data::Dict{String, Any}
    simdir::String
    optimizer                    # stored so fix_investments! can rebuild block models
    ptdf_matrix::Matrix{Float64}
    gen_bus_map::Dict{Int, Int}
    rate_a_nonzero::Set{Int}

    R::Int; N::Int; E::Int; T::Int; T_b::Int; G::Int
    obj_scale::Float64   # multiplied onto op_cost before ADMM penalty; unscaled in get_objective_value

    function PHSubproblem(data::Dict{String, Any}, optimizer, simdir::String;
                                 n_blocks::Int=-1,
                                 max_ph_iterations::Int=-1,
                                 ph_tol::Float64=-1.0,
                                 max_ptdf_iterations::Int=-1,
                                 max_ptdf_per_iteration::Int=-1,
                                 ptdf_tol::Float64=-1.0)
        n_blocks               = n_blocks               == -1   ? get(data["param"], "n_blocks",               2)    : n_blocks
        max_ph_iterations      = max_ph_iterations      == -1   ? get(data["param"], "max_ph_iterations",      100)   : max_ph_iterations
        ph_tol                 = ph_tol                 == -1.0 ? get(data["param"], "ph_tol",                 1e-4) : ph_tol
        max_ptdf_iterations    = max_ptdf_iterations    == -1   ? get(data["param"], "max_ptdf_iterations",    256)  : max_ptdf_iterations
        max_ptdf_per_iteration = max_ptdf_per_iteration == -1   ? get(data["param"], "max_ptdf_per_iteration", 32)   : max_ptdf_per_iteration
        ptdf_tol               = ptdf_tol               == -1.0 ? get(data["param"], "ptdf_tol",               1e-6) : ptdf_tol

        rho       = Float64(get(data["param"], "ph_rho",       1.0))
        obj_scale = Float64(get(data["param"], "ph_obj_scale", 1.0))
        R = data["param"]["num_representatives"]
        N = length(data["bus"])
        E = length(data["branch"])
        T = data["param"]["num_hours"]
        G = length(data["gen"])

        @assert T % n_blocks == 0 "num_hours ($T) must be divisible by n_blocks ($n_blocks)"
        T_b = T ÷ n_blocks

        ptdf_matrix = do_all_ptdf(data)
        if haskey(data["param"], "ptdf_cutoff") && data["param"]["ptdf_cutoff"] != false
            ptdf_matrix = ptdf_matrix .* (map(abs, ptdf_matrix) .>= data["param"]["ptdf_cutoff"])
        end

        _, rate_a_nonzero_vec = get_rate_a_zero(data)
        rate_a_nonzero = Set(parse(Int, x) for x in rate_a_nonzero_vec)
        gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"]
                          for g in keys(data["gen"]))

        block_models = Matrix{PHBlockModel}(undef, R, n_blocks)
        for r in 1:R, b in 1:n_blocks
            t_start = (b - 1) * T_b + 1
            jm = _create_block_jump_model(data, optimizer, r, b, n_blocks,
                                          T_b, t_start, zeros(E), zeros(N))
            block_models[r, b] = PHBlockModel(
                jm, data, simdir,
                ptdf_matrix,
                Dict{Tuple{Int,Int,Int,Bool}, Bool}(),
                rate_a_nonzero,
                max_ptdf_iterations, max_ptdf_per_iteration, ptdf_tol,
                0.0,
                Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}(),
                gen_bus_map,
                1, N, E, T_b, G,  # R=1 per block (single representative)
                b, n_blocks, r, t_start
            )
        end

        n_bd = n_blocks - 1
        lambda  = zeros(N, R, n_bd)
        soc_bar = zeros(N, R, n_bd)
        soc_out = zeros(N, R, n_bd)
        soc_in  = zeros(N, R, n_bd)

        new(block_models, n_blocks,
            nothing, Int[],
            lambda, rho, soc_bar, soc_out, soc_in,
            max_ph_iterations, ph_tol,
            data, simdir, optimizer, ptdf_matrix, gen_bus_map, rate_a_nonzero,
            R, N, E, T, T_b, G, obj_scale)
    end
end

# ---------------------------------------------------------------------------
# _create_block_jump_model
#
# JuMP model for representative r, block b, covering global hours
# t_start : t_start+T_b-1.  gamma and s_energy are injected as fixed
# constants so existing add_constraints! code (which references jump_model[:gamma])
# still works without modification.
# ---------------------------------------------------------------------------
function _create_block_jump_model(data, optimizer, r::Int, block_id::Int,
                                   n_blocks::Int, T_b::Int, t_start::Int,
                                   gamma_val::Vector{Float64},
                                   s_energy_val::Vector{Float64})
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])
    _, rate_a_nonzero = get_rate_a_zero(data)

    η         = data["param"]["bess_efficiency"]
    σ         = get(data["param"], "self_discharge", 0)
    S         = data["param"]["storage_energy_size"]
    soc_ratio = get(data["param"], "soc_init_end_ratio", 0.5)

    m = JuMP.Model(optimizer)

    # --- Generator dispatch ---
    @variable(m, pg[g in 1:G, t in 1:T_b])
    for g in 1:G
        gen = data["gen"]["$g"]
        is_renewable = gen["gen_type"] ∈ data["param"]["renewable_types"]
        is_foreign   = gen["gen_type"] ∈ ["foreign"]
        if is_renewable
            set_lower_bound.(pg[g, :], 0.0)
            for t in 1:T_b
                set_upper_bound(pg[g, t], max(0.0, gen["profile"]["$r"][t_start + t - 1]))
            end
        elseif is_foreign
            for t in 1:T_b
                fix(pg[g, t], gen["profile"]["$r"][t_start + t - 1]; force=true)
            end
        else
            set_lower_bound.(pg[g, :], gen["pmin"])
            set_upper_bound.(pg[g, :], gen["pmax"])
        end
    end

    # gamma fixed as a JuMP variable so add_constraints! works without changes
    @variable(m, gamma[a=1:E])
    for a in 1:E
        fix(m[:gamma][a], gamma_val[a]; force=true)
    end

    # --- Storage variables ---
    @variable(m, ue[i=1:N, t=1:T_b] >= 0)
    @variable(m, soc[i=1:N, t=1:T_b] >= 0)
    @variable(m, ch[i=1:N, t=1:T_b] >= 0)
    @variable(m, dis[i=1:N, t=1:T_b] >= 0)

    for i in 1:N
        cap = s_energy_val[i] * S
        set_upper_bound.(soc[i, :], max(cap, 0.0))
        set_upper_bound.(ch[i, :],  max(cap / 4, 0.0))
        set_upper_bound.(dis[i, :], max(cap / 4, 0.0))
    end

    # --- Constraints ---
    gen_bus_map = Dict(parse(Int, g) => data["gen"]["$g"]["gen_bus"] for g in keys(data["gen"]))

    @constraint(m,
        power_balance[t in 1:T_b],
        sum(pg[g, t] for g in 1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t_start + t - 1] for i in 1:N)
        + sum(ue[i, t] for i in 1:N)
        - sum(ch[i, t] for i in 1:N)
        + sum(dis[i, t] for i in 1:N)
        == 0
    )

    @constraint(m,
        soc_over_time[i in 1:N, t in 2:T_b],
        soc[i, t] == (1 - σ) * soc[i, t-1] + ch[i, t] * η - dis[i, t] / η
    )

    @constraint(m, soc_energy_ub[i in 1:N, t in 1:T_b], soc[i, t] <= s_energy_val[i] * S)
    @constraint(m, charge_ub[i in 1:N, t in 1:T_b],     ch[i, t]  <= s_energy_val[i] * S / 4)
    @constraint(m, discharge_ub[i in 1:N, t in 1:T_b],  dis[i, t] <= s_energy_val[i] * S / 4)

    # --- SOC boundary conditions ---
    if block_id == 1
        # Left boundary: fixed initial SOC from soc_ratio parameter
        @constraint(m,
            soc_start[i in 1:N],
            soc[i, 1] == soc_ratio * s_energy_val[i] * S + ch[i, 1] * η - dis[i, 1] / η
        )
    else
        # Left boundary: soc_in_b[i] is the ADMM linking variable from the previous block
        @variable(m, soc_in_b[i=1:N] >= 0)
        for i in 1:N
            set_upper_bound(m[:soc_in_b][i], max(s_energy_val[i] * S, 0.0))
        end
        @constraint(m,
            soc_start[i in 1:N],
            soc[i, 1] == (1 - σ) * m[:soc_in_b][i] + ch[i, 1] * η - dis[i, 1] / η
        )
    end

    if block_id == n_blocks
        # Right boundary: terminal SOC matches soc_ratio parameter
        @constraint(m,
            soc_end[i in 1:N],
            soc[i, T_b] == soc_ratio * s_energy_val[i] * S
        )
    end
    # Interior and first blocks have no soc_end — soc[i, T_b] is the free right linking variable

    # Container for lazy PTDF constraints
    m[:ptdf_flow] = Dict{String, ConstraintRef}()

    # Placeholder objective — replaced by set_block_objective! before each solve
    @objective(m, Min, 0.0)

    return m
end

# ---------------------------------------------------------------------------
# PHHub
#
# Persistent object (one per worker in the distributed framework) that
# tracks which PTDF lazy constraints were active across Benders iterations.
# The PHSubproblem itself is rebuilt each Benders iteration (because
# investments change), but the Hub carries forward the discovered constraint
# index sets so block models can be warm-started instead of rediscovering
# them from scratch.
#
# Only constraint indices (a, t, ub) are stored — not JuMP ConstraintRef
# objects — so this is fully serialization-safe and investment-agnostic.
# The actual constraint expressions (which depend on gamma_val) are rebuilt
# from scratch during warm_start_ph! using the current gamma values.
# ---------------------------------------------------------------------------
mutable struct PHHub
    # tracked_constraints[r, b]: set of (a, r_in_block, t, ub) keys that were
    # active at the end of the last PH solve for rep r, block b.
    tracked_constraints::Matrix{Dict{Tuple{Int,Int,Int,Bool}, Bool}}

    # constraint_stats[r, b]: violation history — useful for pruning stale
    # constraints across Benders iterations.
    constraint_stats::Matrix{Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}}

    n_blocks::Int
    R::Int
    T_b::Int
    T::Int

    function PHHub(R::Int, n_blocks::Int, T_b::Int, T::Int)
        tc = [Dict{Tuple{Int,Int,Int,Bool}, Bool}() for _ in 1:R, _ in 1:n_blocks]
        cs = [Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}() for _ in 1:R, _ in 1:n_blocks]
        new(tc, cs, n_blocks, R, T_b, T)
    end
end

# ---------------------------------------------------------------------------
# fix_investments!
#
# Rebuilds all block models with the given master investment values and
# resets all ADMM state to zero.  Pass hub= to warm-start PTDF constraints
# from the previous Benders iteration.
# ---------------------------------------------------------------------------
function fix_investments!(ph::PHSubproblem, gamma_val::Vector, s_energy_val::Vector;
                          hub::Union{PHHub, Nothing}=nothing)
    for r in 1:ph.R, b in 1:ph.n_blocks
        t_start = (b - 1) * ph.T_b + 1
        jm = _create_block_jump_model(ph.data, ph.optimizer, r, b, ph.n_blocks,
                                       ph.T_b, t_start, gamma_val, s_energy_val)
        bm = ph.block_models[r, b]
        bm.jump_model          = jm
        bm.tracked_constraints = Dict{Tuple{Int,Int,Int,Bool}, Bool}()
        bm.constraint_stats    = Dict{Tuple{Int,Int,Int,Bool}, ConstraintStat}()
        bm.solve_time          = 0.0
    end

    if !isnothing(hub)
        warm_start_ph!(ph, hub, gamma_val)
    end

    fill!(ph.lambda,  0.0)
    fill!(ph.soc_bar, 0.0)
    fill!(ph.soc_out, 0.0)
    fill!(ph.soc_in,  0.0)

    ph.storage_nodes       = findall(s_energy_val .> 0)
    ph.current_investments = (copy(gamma_val), copy(s_energy_val))
    return nothing
end

# ---------------------------------------------------------------------------
# set_block_objective!
#
# Sets the ADMM-augmented quadratic objective for block b of representative r.
#
# Each block contributes penalties from its adjacent boundaries:
#   Left boundary  (bd = b-1, if b > 1):
#       - λ[i,r,bd] * soc_in_b[i]  +  (ρ/2) * (soc_out[i,r,bd] - soc_in_b[i])²
#   Right boundary (bd = b,   if b < n_blocks):
#       + λ[i,r,bd] * soc[i,T_b]   +  (ρ/2) * (soc[i,T_b] - soc_bar[i,r,bd])²
#
# Produces a QP — Gurobi handles this natively.
# ---------------------------------------------------------------------------
function set_block_objective!(ph::PHSubproblem, r::Int, b::Int)
    bm   = ph.block_models[r, b]
    m    = bm.jump_model
    data = ph.data
    N, T_b, G = ph.N, ph.T_b, ph.G

    prob_r = data["param"]["representative_prob"][r]
    op_wt  = get(data["param"], "operational_weight", 1)
    stor_op_cost = get(data["param"], "storage_operation_cost", 0.0)

    op_cost = ph.obj_scale * prob_r * op_wt * (
        sum(compute_gen_cost(m[:pg][g, t], data["gen"]["$g"]) for g in 1:G, t in 1:T_b) +
        sum(data["param"]["under_served_penalty"] * m[:ue][i, t] for i in 1:N, t in 1:T_b) +
        sum(stor_op_cost * (m[:ch][i, t] + m[:dis][i, t]) for i in 1:N, t in 1:T_b)
    )

    # QuadExpr accumulates both linear (λ) and quadratic (ρ) penalty terms.
    # With rho=0 this stays an LP — useful for validating correctness first.
    penalty = zero(QuadExpr)

    # Left boundary penalty (block b receives SOC from block b-1)
    if b > 1
        bd = b - 1
        soc_in_b = m[:soc_in_b]
        for i in 1:N
            add_to_expression!(penalty,
                -ph.lambda[i, r, bd], soc_in_b[i]
            )
            if ph.rho > 0
                # (soc_out_fixed - soc_in_b)^2 = soc_in_b^2 - 2*soc_out_fixed*soc_in_b + const
                add_to_expression!(penalty, ph.rho / 2, soc_in_b[i], soc_in_b[i])
                add_to_expression!(penalty, -ph.rho * ph.soc_out[i, r, bd], soc_in_b[i])
            end
        end
    end

    # Right boundary penalty (block b hands SOC to block b+1)
    if b < ph.n_blocks
        bd = b
        for i in 1:N
            add_to_expression!(penalty,
                ph.lambda[i, r, bd], m[:soc][i, T_b]
            )
            if ph.rho > 0
                # (soc_out - soc_bar)^2 = soc_out^2 - 2*soc_bar*soc_out + const
                add_to_expression!(penalty, ph.rho / 2, m[:soc][i, T_b], m[:soc][i, T_b])
                add_to_expression!(penalty, -ph.rho * ph.soc_bar[i, r, bd], m[:soc][i, T_b])
            end
        end
    end

    @objective(m, Min, op_cost + penalty)
    return nothing
end

# ---------------------------------------------------------------------------
# PTDF helpers — multiple dispatch on PHBlockModel
#
# Variables in PHBlockModel are indexed [node/gen, t] (no leading r dimension
# since each PHBlockModel covers exactly one representative). Load profiles are
# offset by t_start so they refer to the correct global hours.
# ---------------------------------------------------------------------------
function find_violations(bm::PHBlockModel, flows::Matrix{Float64},
                          gamma_values::Vector{Float64}, niter::Int)
    violations = []
    for a in bm.rate_a_nonzero
        cap_increment    = get_capacity_increment(bm.data, a)
        line_limit_fixed = bm.data["branch"]["$a"]["rate_a"] + gamma_values[a] * cap_increment

        for t in 1:bm.T
            flow_val         = flows[a, t]
            violation_amount = abs(flow_val) - line_limit_fixed
            ub  = (flow_val >= 0)
            key = (a, 1, t, ub)  # r dimension collapsed to 1 within a block model

            if haskey(bm.constraint_stats, key)
                stat = bm.constraint_stats[key]
                push!(stat.tightness_history, (niter, violation_amount))
                if violation_amount > bm.ptdf_tol
                    stat.total_appearances += 1
                    stat.max_violation = max(stat.max_violation, violation_amount)
                end
            elseif violation_amount > bm.ptdf_tol
                bm.constraint_stats[key] = ConstraintStat(
                    niter, 1, violation_amount, [(niter, violation_amount)]
                )
            end

            get(bm.tracked_constraints, key, false) && continue
            violation_amount > bm.ptdf_tol && push!(violations, (a, 1, t, ub, violation_amount))
        end
    end
    return violations
end

function add_constraints!(bm::PHBlockModel, sorted_violations, gamma_val::Vector{Float64})
    n_added = 0
    for (a, _, t, ub, _) in Iterators.take(sorted_violations, bm.max_ptdf_per_iteration)
        ptdf_row        = bm.ptdf_matrix[a, :]
        cap_increment   = get_capacity_increment(bm.data, a)
        line_limit_expr = bm.data["branch"]["$a"]["rate_a"] + bm.jump_model[:gamma][a] * cap_increment

        flow_expr = sum(ptdf_row[i] * (
            sum(bm.jump_model[:pg][g, t]
                for g in 1:bm.G if bm.gen_bus_map[g] == i; init=0.0)
            - bm.data["bus"]["$i"]["load"]["$(bm.r)"][bm.t_start + t - 1]
            + bm.jump_model[:ue][i, t]
            - bm.jump_model[:ch][i, t]
            + bm.jump_model[:dis][i, t]
        ) for i in 1:bm.N)

        ckey = ub ? "$(a)_$(t)_ub" : "$(a)_$(t)_lb"
        if ub
            bm.jump_model[:ptdf_flow][ckey] = @constraint(bm.jump_model, flow_expr <= line_limit_expr)
        else
            bm.jump_model[:ptdf_flow][ckey] = @constraint(bm.jump_model, flow_expr >= -line_limit_expr)
        end
        bm.tracked_constraints[(a, 1, t, ub)] = true
        n_added += 1
    end
    return n_added
end

function prune_stale_constraints!(bm::PHBlockModel, stale_k::Int)
    to_prune = Tuple{Int,Int,Int,Bool}[]

    for (key, _) in bm.tracked_constraints
        stat = get(bm.constraint_stats, key, nothing)
        isnothing(stat) && continue
        history = stat.tightness_history
        length(history) < stale_k && continue
        if stat.total_appearances == 1 && all(v <= bm.ptdf_tol for (_, v) in history[end-stale_k+1:end])
            push!(to_prune, key)
        end
    end

    for key in to_prune
        (a, _, t, ub) = key
        ckey = ub ? "$(a)_$(t)_ub" : "$(a)_$(t)_lb"
        if haskey(bm.jump_model[:ptdf_flow], ckey)
            delete(bm.jump_model, bm.jump_model[:ptdf_flow][ckey])
            delete!(bm.jump_model[:ptdf_flow], ckey)
        end
        delete!(bm.tracked_constraints, key)
    end

    return length(to_prune)
end

function save_tracked_constraints(bm::PHBlockModel, n_violated::Int)
    ks = keys(bm.tracked_constraints)
    df = DataFrame(
        arc             = [k[1] for k in ks],
        time            = [k[3] for k in ks],
        ub              = [k[4] for k in ks],
        tracked         = collect(values(bm.tracked_constraints)),
        violations_left = fill(n_violated, length(ks)),
    )
    output_dir = joinpath(bm.simdir, "output")
    mkpath(output_dir)
    CSV.write(joinpath(output_dir, "ph_tracked_constraints_r$(bm.r)_b$(bm.block_id).csv"), df)
end

function compute_block_flows!(flows::Matrix{Float64}, bm::PHBlockModel)
    pg_values  = value.(bm.jump_model[:pg])
    ue_values  = value.(bm.jump_model[:ue])
    ch_values  = value.(bm.jump_model[:ch])
    dis_values = value.(bm.jump_model[:dis])
    r_str      = "$(bm.r)"

    net_injections = zeros(bm.N)
    for t in 1:bm.T
        fill!(net_injections, 0.0)
        for i in 1:bm.N
            net_injections[i] = sum(pg_values[g, t]
                for g in 1:bm.G if bm.gen_bus_map[g] == i; init=0.0)
            net_injections[i] -= bm.data["bus"]["$i"]["load"][r_str][bm.t_start + t - 1]
            net_injections[i] += ue_values[i, t]
            net_injections[i] -= ch_values[i, t]
            net_injections[i] += dis_values[i, t]
        end
        flows[:, t] = bm.ptdf_matrix * net_injections
    end
end

# PTDF lazy-constraint loop for a single PHBlockModel
function solve_block_ptdf!(bm::PHBlockModel, gamma_val::Vector{Float64})
    solved = false
    niter  = 0
    t0     = time()
    flows  = zeros(bm.E, bm.T)

    while !solved && niter < bm.max_ptdf_iterations
        JuMP.optimize!(bm.jump_model)
        st = termination_status(bm.jump_model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || (println("  Block $(bm.block_id) rep $(bm.r): solver status $st"); break)

        compute_block_flows!(flows, bm)
        violations        = find_violations(bm, flows, gamma_val, niter + 1)
        sorted_violations = sort(violations, by = x -> -x[5])
        n_added  = add_constraints!(bm, sorted_violations, gamma_val)
        stale_k  = get(bm.data["param"], "stale_k", nothing)
        n_pruned = (isnothing(stale_k) || isempty(sorted_violations)) ? 0 : prune_stale_constraints!(bm, stale_k)

        if n_added > 0 || n_pruned > 0
            save_tracked_constraints(bm, length(sorted_violations))
        end

        solved = isempty(sorted_violations)
        niter += 1
        println("  Block $(bm.block_id) rep $(bm.r) | PTDF iter $niter: $(length(sorted_violations)) violations, $n_added added, $n_pruned pruned")
    end

    bm.solve_time += time() - t0
    return termination_status(bm.jump_model)
end

# ---------------------------------------------------------------------------
# solve!
#
# Outer ADMM loop. Each iteration:
#   1. Set augmented objectives for all blocks.
#   2. Solve blocks sequentially b=1:n_blocks (Gauss-Seidel order — each
#      block sees the freshest boundary values from the block just before it).
#   3. After the full pass, update consensus and multipliers.
#   4. Check convergence on primal residual ||soc_out - soc_in||_inf.
#
# Representatives are independent within each block solve and can be
# parallelised with Threads.@threads when ready.
# ---------------------------------------------------------------------------
function solve!(ph::PHSubproblem)
    isnothing(ph.current_investments) && error("Call fix_investments! before solve!")
    gamma_val, _ = ph.current_investments

    soc_log_rows = NamedTuple{(:iteration, :rep, :boundary, :node, :soc_out, :soc_in),
                               Tuple{Int,Int,Int,Int,Float64,Float64}}[]
    conv_log_rows = NamedTuple{(:iteration, :iter_time_s, :elapsed_time_s, :mean_residual, :max_residual),
                                Tuple{Int,Float64,Float64,Float64,Float64}}[]

    t_solve_start = time()

    for ph_iter in 1:ph.max_ph_iterations
        println("=== PH/ADMM iteration $ph_iter ===")
        t_iter_start = time()

        for b in 1:ph.n_blocks
            for r in 1:ph.R
                set_block_objective!(ph, r, b)
                solve_block_ptdf!(ph.block_models[r, b], gamma_val)
            end

            # Extract boundary values immediately after each block so the next
            # block sees the freshest values (Gauss-Seidel)
            for r in 1:ph.R
                bm = ph.block_models[r, b]
                # Right boundary of block b  →  soc_out[:, r, b]
                if b < ph.n_blocks
                    for i in 1:ph.N
                        ph.soc_out[i, r, b] = value(bm.jump_model[:soc][i, ph.T_b])
                    end
                end
                # Left boundary of block b   →  soc_in[:, r, b-1]
                if b > 1
                    for i in 1:ph.N
                        ph.soc_in[i, r, b-1] = value(bm.jump_model[:soc_in_b][i])
                    end
                end
            end
        end

        update_admm!(ph)

        # Log boundary SOC for storage nodes after each ADMM iteration
        if ph.n_blocks > 1 && !isempty(ph.storage_nodes)
            for bd in 1:(ph.n_blocks - 1), r in 1:ph.R, i in ph.storage_nodes
                push!(soc_log_rows, (iteration=ph_iter, rep=r, boundary=bd, node=i,
                                     soc_out=ph.soc_out[i, r, bd], soc_in=ph.soc_in[i, r, bd]))
            end
            CSV.write(joinpath(ph.simdir, "output", "soc_boundary_log.csv"), DataFrame(soc_log_rows))
        end

        residuals = ph.n_blocks > 1 ? abs.(ph.soc_out .- ph.soc_in) : zeros(0)
        max_residual  = isempty(residuals) ? 0.0 : maximum(residuals)
        mean_residual = isempty(residuals) ? 0.0 : sum(residuals) / length(residuals)
        iter_time     = time() - t_iter_start
        elapsed_time  = time() - t_solve_start

        push!(conv_log_rows, (iteration=ph_iter, iter_time_s=iter_time,
                               elapsed_time_s=elapsed_time,
                               mean_residual=mean_residual, max_residual=max_residual))
        CSV.write(joinpath(ph.simdir, "output", "admm_convergence_log.csv"), DataFrame(conv_log_rows))

        println("PH/ADMM iter $ph_iter | max residual = $max_residual | mean residual = $mean_residual")

        max_residual < ph.ph_tol && (println("Converged in $ph_iter iterations."); break)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# update_admm!
#
# Consensus:   soc_bar = (soc_out + soc_in) / 2
# Multipliers: λ      += ρ  * (soc_out - soc_in)
#
# Both broadcasts operate over [N, R, n_blocks-1] simultaneously.
# ---------------------------------------------------------------------------
function update_admm!(ph::PHSubproblem)
    @. ph.soc_bar = (ph.soc_out + ph.soc_in) / 2
    @. ph.lambda  = ph.lambda + ph.rho * (ph.soc_out - ph.soc_in)
    return nothing
end

function get_objective_value(ph::PHSubproblem)
    # objective_value includes the ADMM penalty and obj_scale — strip both to get true op cost
    # instead, recompute from primal values directly
    data = ph.data
    op_wt        = get(data["param"], "operational_weight", 1)
    stor_op_cost = get(data["param"], "storage_operation_cost", 0.0)
    total = 0.0
    for r in 1:ph.R, b in 1:ph.n_blocks
        bm     = ph.block_models[r, b]
        m      = bm.jump_model
        prob_r = data["param"]["representative_prob"][r]
        pg_v   = value.(m[:pg])
        ue_v   = value.(m[:ue])
        ch_v   = value.(m[:ch])
        dis_v  = value.(m[:dis])
        for g in 1:ph.G, t in 1:bm.T
            total += prob_r * op_wt * compute_gen_cost(pg_v[g, t], data["gen"]["$g"])
        end
        for i in 1:ph.N, t in 1:bm.T
            total += prob_r * op_wt * data["param"]["under_served_penalty"] * ue_v[i, t]
            total += prob_r * op_wt * stor_op_cost * (ch_v[i, t] + dis_v[i, t])
        end
    end
    return total
end

# Convenience constructor from a PHSubproblem (after first solve)
function PHHub(ph::PHSubproblem)
    hub = PHHub(ph.R, ph.n_blocks, ph.T_b, ph.T)
    sync_from_ph!(hub, ph)
    return hub
end

# ---------------------------------------------------------------------------
# sync_from_ph!
#
# Called after each PH solve completes. Copies the tracked constraint indices
# and violation stats from all block models into the hub for the next
# Benders iteration's warm-start.
# ---------------------------------------------------------------------------
function sync_from_ph!(hub::PHHub, ph::PHSubproblem)
    for r in 1:ph.R, b in 1:ph.n_blocks
        bm = ph.block_models[r, b]
        hub.tracked_constraints[r, b] = copy(bm.tracked_constraints)
        hub.constraint_stats[r, b]    = copy(bm.constraint_stats)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# warm_start_ph!
#
# Called inside fix_investments! (when a hub is provided) after block models
# are freshly built. Re-injects all previously tracked constraints into each
# block model using the current gamma_val for line limits.
#
# Constraints that are slack under the new investments will simply not bind —
# this is safe, just slightly larger models. Constraints that are no longer
# valid (gamma reduced below tracked value) are still added conservatively;
# the PTDF loop will verify feasibility and prune if needed.
# ---------------------------------------------------------------------------
function warm_start_ph!(ph::PHSubproblem, hub::PHHub, gamma_val::Vector{Float64})
    total_added = 0
    for r in 1:ph.R, b in 1:ph.n_blocks
        bm = ph.block_models[r, b]
        n_added = _inject_tracked_constraints!(bm, hub.tracked_constraints[r, b], gamma_val)
        # Carry forward stats so tightness history is preserved
        bm.constraint_stats = copy(hub.constraint_stats[r, b])
        total_added += n_added
    end
    println("PHHub warm-start: injected $total_added constraints across $(ph.R * ph.n_blocks) block models")
    return nothing
end

# ---------------------------------------------------------------------------
# _inject_tracked_constraints!
#
# Adds a set of constraint indices into a block model. Builds the flow
# expression from scratch using the current gamma_val — this is what makes
# warm-starting investment-agnostic. Skips any keys already present.
# ---------------------------------------------------------------------------
function _inject_tracked_constraints!(bm::PHBlockModel,
                                       keys::Dict{Tuple{Int,Int,Int,Bool}, Bool},
                                       gamma_val::Vector{Float64})
    n_added = 0
    for (key, _) in keys
        (a, _, t, ub) = key
        get(bm.tracked_constraints, key, false) && continue  # already present

        ptdf_row        = bm.ptdf_matrix[a, :]
        cap_increment   = get_capacity_increment(bm.data, a)
        line_limit_expr = bm.data["branch"]["$a"]["rate_a"] + bm.jump_model[:gamma][a] * cap_increment

        flow_expr = sum(ptdf_row[i] * (
            sum(bm.jump_model[:pg][g, t]
                for g in 1:bm.G if bm.gen_bus_map[g] == i; init=0.0)
            - bm.data["bus"]["$i"]["load"]["$(bm.r)"][bm.t_start + t - 1]
            + bm.jump_model[:ue][i, t]
            - bm.jump_model[:ch][i, t]
            + bm.jump_model[:dis][i, t]
        ) for i in 1:bm.N)

        ckey = ub ? "$(a)_$(t)_ub" : "$(a)_$(t)_lb"
        if ub
            bm.jump_model[:ptdf_flow][ckey] = @constraint(bm.jump_model, flow_expr <= line_limit_expr)
        else
            bm.jump_model[:ptdf_flow][ckey] = @constraint(bm.jump_model, flow_expr >= -line_limit_expr)
        end
        bm.tracked_constraints[key] = true
        n_added += 1
    end
    return n_added
end