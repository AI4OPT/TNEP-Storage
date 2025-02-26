using JuMP
using CSV, DataFrames
include("../../../helpers/compute_gen_cost.jl")
include("../storage_candidates/naive_candidates.jl")
include("../rate_a_zero.jl")
include("base_ptdf.jl")

function create_model_r1_ptdf_iterative(simdir, data::Dict{String, Any}, optimizer; 
    line_investments=nothing, 
    storage_investments=nothing,
    max_ptdf_iterations::Int=128,
    max_ptdf_per_iteration::Int=64,
    ptdf_tol::Float64=1e-6)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]

    # Initialize model and basic components
    model = create_base_model(data, optimizer, line_investments=line_investments, storage_investments=storage_investments)
    # Set up model extensions for metadata
    model.ext[:solve_metadata] = Dict(
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
        :ptdf_tol => ptdf_tol
    )
    logfile = joinpath(simdir, "gurobi_logfile.log")
    set_optimizer_attribute(model, "LogFile", logfile)
    set_optimizer_attribute(model, "MIPGap", data["param"]["mip_gap"])

    # Store PTDF matrix in model extension
    model.ext[:PTDF] = do_all_ptdf(data)
    
    # Track which branches already have constraints
    model.ext[:tracked_constraints] = Dict{Tuple{Int,Int,Int}, Bool}()

    # Initialize rate_a lookup
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)
    model.ext[:rate_a_nonzero] = Set(parse(Int, x) for x in rate_a_nonzero)
    
    # Add container for flow constraints
    model[:ptdf_flow] = Dict{String, ConstraintRef}()

    # warm-start the PTDF with certain flow constraints
    # Add after PTDF matrix storage
    if isfile(joinpath(simdir, "congested_constraints.csv"))
        constraints_df = CSV.read(joinpath(simdir, "congested_constraints.csv"), DataFrame)
        
        rows_done = 0
        for row in eachrow(constraints_df)
            rows_done += 1
            println(rows_done)
            a, r, t = row.arc, row.rep, row.time
            if row.tracked
                line_limit = data["branch"]["$a"]["rate_a"] + model[:gamma][a] * get_capacity_increment(data, a)
                row = model.ext[:PTDF][a,:]
                
                model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(model,
                    sum(row[i] * (
                        sum(model[:pg][r,g,t] for g in 1:G if data["gen"]["$g"]["gen_bus"] == i; init=0.0)
                        - data["bus"]["$i"]["load"]["$r"][t]
                        + model[:ue][r,i,t] - model[:ch][r,i,t]
                    ) for i in 1:length(data["bus"])) <= line_limit
                )
                
                model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(model,
                    sum(row[i] * (
                        sum(model[:pg][r,g,t] for g in 1:G if data["gen"]["$g"]["gen_bus"] == i; init=0.0)
                        - data["bus"]["$i"]["load"]["$r"][t]
                        + model[:ue][r,i,t] - model[:ch][r,i,t]
                    ) for i in 1:length(data["bus"])) >= -line_limit
                )
                
                model.ext[:tracked_constraints][(a,r,t)] = true
            end
        end
    end

    # Begin lazy PTDF loop
    solved = false
    niter = 0
    t0 = time()

    while !solved && niter < model.ext[:solve_metadata][:max_ptdf_iterations]
        # Solve model
        optimize!(model)
        
        # Exit if not solved optimally
        st = termination_status(model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get current solution values
        pg_values = value.(model[:pg])
        ue_values = value.(model[:ue])
        ch_values = value.(model[:ch])
        gamma_values = value.(model[:gamma])
        
        # Compute all flows at once
        flows = zeros(length(model.ext[:rate_a_nonzero]), R, T)
        compute_flows!(flows, pg_values, ue_values, ch_values, data, model.ext[:PTDF])
        
        # Check for violations
        n_violated = 0
        n_added = 0
        
        # First find violations without modifying the model
        checked_count = 0
        for a in model.ext[:rate_a_nonzero]
            checked_count += 1
            line_limit = data["branch"]["$a"]["rate_a"] + gamma_values[a] * get_capacity_increment(data, a)
            
            for r in 1:R, t in 1:T
                if get(model.ext[:tracked_constraints], (a,r,t), false)
                    continue
                end
                
                if abs(flows[a,r,t]) > line_limit + model.ext[:solve_metadata][:ptdf_tol]
                    n_violated += 1
                    
                    if n_added < model.ext[:solve_metadata][:max_ptdf_per_iteration]
                        row = model.ext[:PTDF][a,:]
                        
                        # Add both constraints
                        model[:ptdf_flow]["$(a)_$(r)_$(t)_ub"] = @constraint(model,
                            sum(row[i] * (
                                sum(model[:pg][r,g,t] for g in 1:G if data["gen"]["$g"]["gen_bus"] == i; init=0.0)
                                - data["bus"]["$i"]["load"]["$r"][t]
                                + model[:ue][r,i,t] - model[:ch][r,i,t]
                            ) for i in 1:length(data["bus"])) <= data["branch"]["$a"]["rate_a"] + model[:gamma][a] * get_capacity_increment(data, a)
                        )
                        
                        model[:ptdf_flow]["$(a)_$(r)_$(t)_lb"] = @constraint(model,
                            sum(row[i] * (
                                sum(model[:pg][r,g,t] for g in 1:G if data["gen"]["$g"]["gen_bus"] == i; init=0.0)
                                - data["bus"]["$i"]["load"]["$r"][t]
                                + model[:ue][r,i,t] - model[:ch][r,i,t]
                            ) for i in 1:length(data["bus"])) >= -(data["branch"]["$a"]["rate_a"] + model[:gamma][a] * get_capacity_increment(data, a))
                        )
                        
                        model.ext[:tracked_constraints][(a,r,t)] = true
                        n_added += 1
                        
                        
                        if n_added >= model.ext[:solve_metadata][:max_ptdf_per_iteration]
                            break
                        end
                    end
                end
            end
            
            if n_added >= model.ext[:solve_metadata][:max_ptdf_per_iteration]
                break
            end
        end
        
        solved = (n_violated == 0)
        niter += 1
        println("approx. progress $(checked_count / E)")
        println("Iteration $niter: Found $n_violated violations, added $n_added constraints")
    end

    # Record solve time
    solve_time = time() - t0
    model.ext[:solve_metadata][:solve_time] = solve_time

    save_tracked_constraints(simdir, model)
    save_power_injections(simdir, model, data)
    return model
end

function create_base_model(data::Dict{String, Any}, optimizer; 
    line_investments=nothing, 
    storage_investments=nothing)

    # Initialize sets
    R = data["param"]["num_representatives"]
    N = length(data["bus"])
    E = length(data["branch"])
    T = data["param"]["num_hours"]
    G = length(data["gen"])
    K = data["param"]["num_cap_upgrades_max"]
    rate_a_zero, rate_a_nonzero = get_rate_a_zero(data)

    model = JuMP.Model(optimizer)

    #
    #   I. Variables
    #
    
    # generator active dispatch
    nonrenewable_generators = filter(g -> lowercase(data["gen"][g]["gen_type"]) ∉ data["param"]["renewable_types"], keys(data["gen"]))
    JuMP.@variable(model, pg[r in 1:R, g in 1:G, t in 1:T])

    for g in 1:G
        gen = data["gen"]["$g"]
        is_renewable = gen["gen_type"] ∈ data["param"]["renewable_types"]
        is_foreign = gen["gen_type"]  ∈ ["foreign"]
        if is_renewable
            set_lower_bound.(pg[:, g, :], 0.0)
            for r in 1:R, t in 1:T
                set_upper_bound(pg[r, g, t], max(0, gen["profile"]["$r"][t]))
            end
        elseif is_foreign
            for r in 1:R, t in 1:T
                fix(pg[r, g, t], gen["profile"]["$r"][t])
            end
        else
            set_lower_bound.(pg[:, g, :], gen["pmin"])
            set_upper_bound.(pg[:, g, :], gen["pmax"])
        end
    end

    JuMP.@variable(model, ue[r=1:R, i=1:N, t=1:T] >= 0) # under-served energy at bus
    JuMP.@variable(model, 0 <= gamma[a=1:E] <= K, Int) # investment level of capacity upgrade
    # JuMP.@variable(model, pf[r=1:R, a=1:E, t=1:T]) # branch flows
    JuMP.@variable(model, s_power[i=1:N] >= 0) # power rating of storage
    JuMP.@variable(model, s_energy[i=1:N] >= 0) # energy rating of storage
    JuMP.@variable(model, soc[r=1:R, i=1:N, t=1:T] >= 0) # state of charge of storage
    JuMP.@variable(model, ch[r=1:R, i=1:N, t=1:T]) # charging of storage
    JuMP.@variable(model, sigma[i=1:N], Bin) # binary variable for installation of storage

    # Add all non-PTDF constraints
    # Global power balance
    @constraint(model, 
        power_balance[r in 1:R, t in 1:T],
        sum(pg[r,g,t] for g in 1:G)
        - sum(data["bus"]["$i"]["load"]["$r"][t] for i in 1:N)
        + sum(ue[r,i,t] for i in 1:N) 
        - sum(ch[r,i,t] for i in 1:N) 
        == 0
    )

    # Previous investment constraints if applicable
    nonzero_storage_nodes = Set()
    if haskey(data["param"], "prev_simdir")
        prev_simdir = data["param"]["prev_simdir"]
        line_inv = CSV.read(joinpath(prev_simdir, "output", "line_investments.csv"), DataFrame)
        storage_inv = CSV.read(joinpath(prev_simdir, "output", "storage_investments.csv"), DataFrame)
        for i in 1:N
            if storage_inv[i, :Storage_Power] != 0 || storage_inv[i, :Storage_Energy] != 0
                push!(nonzero_storage_nodes, "$i")
            end
        end
        
        @constraint(model,
            old_gamma[a in 1:E],
            gamma[a] >= line_inv[a, :Upgrade_Lvl]
        )
        @constraint(model,
            old_s_power[i in 1:N],
            s_power[i] >= storage_inv[i, :Storage_Power]
        )
        @constraint(model,
            old_s_energy[i in 1:N],
            s_energy[i] >= storage_inv[i, :Storage_Energy]
        )
    end

    # Investment test constraints if applicable
    if line_investments !== nothing
        inv = CSV.read(line_investments, DataFrame)
        @constraint(model,
            inv_gamma[a in 1:E],
            gamma[a] == inv[a, :Upgrade_Lvl]
        )
    end
    if storage_investments !== nothing
        inv = CSV.read(storage_investments, DataFrame)
        @constraint(model,
            inv_s_power[i in 1:N],
            s_power[i] == inv[i, :Storage_Power]
        )
        @constraint(model,
            inv_s_energy[i in 1:N],
            s_energy[i] == inv[i, :Storage_Energy]
        )
        @constraint(model,
            inv_sigma[i in 1:N],
            sigma[i] == ((inv[i, :Storage_Power] != 0 || inv[i, :Storage_Energy] != 0) ? 1 : 0)
        )
    end

    # if rate a is zero (unlimited), then don't allow upgrades
    JuMP.@constraint(model, 
        rate_a_zero_line_upgrade[a in rate_a_zero],
        gamma[a] == 0
    )

    # soc over time constraint
    JuMP.@constraint(model, 
        soc_over_time[r in 1:R, i in 1:N, t in 2:T],
        soc[r,i,t] == soc[r,i,t-1] + ch[r,i,t]
    )

    # OPTIONAL: soc 0.5 constraint
    JuMP.@constraint(model,
        soc_start[r in 1:R, i in 1:N],
        soc[r,i,1] == 0.5 * s_energy[i] + ch[r,i,1]
    )
    JuMP.@constraint(model,
        soc_end[r in 1:R, i in 1:N],
        soc[r,i,T] == 0.5 * s_energy[i]
    )

    # soc energy rating constraint
    JuMP.@constraint(model, 
        soc_energy_ub[r in 1:R, i in 1:N, t in 1:T],
        soc[r,i,t] <= s_energy[i]
    )

    # energy rating only if storage installed
    JuMP.@constraint(model, 
        installed_energy_ub[i in 1:N],
        s_energy[i] <= sigma[i] * data["param"]["max_energy_rating"]
    )

    # power rating only if storage installed
    JuMP.@constraint(model, 
        installed_power_ub[i in 1:N],
        s_power[i] <= sigma[i] * data["param"]["max_power_rating"]
    )

    # ensure that all storage is short-duration, i.e. can only store 4-hours worth of discharge
    if storage_investments == nothing 
        JuMP.@constraint(model, 
            short_duration[i in 1:N],
            s_energy[i] <= 4.0 * s_power[i]
        )
    end

    # charge/discharge must be constrained by power rating
    JuMP.@constraint(model,
        charge_discharge_lb[r in 1:R, i in 1:N, t in 1:T],
        -s_power[i] <= ch[r,i,t]
    )
    JuMP.@constraint(model,
        charge_discharge_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= s_power[i]
    )

    # charge/discharge must be constrained by installation
    JuMP.@constraint(model,
        charge_discharge_i_lb[r in 1:R, i in 1:N, t in 1:T],
        -sigma[i] * data["param"]["max_power_rating"] <= ch[r,i,t]
    )
    JuMP.@constraint(model,
        charge_discharge_i_ub[r in 1:R, i in 1:N, t in 1:T],
        ch[r,i,t] <= sigma[i] * data["param"]["max_power_rating"]
    )

    # OPTIONAL: CANDIDATE STORAGE LOCATIONS ONLY
    if haskey(data["param"], "candidate_no_upgrades_dir")
        no_upgrades_dir = data["param"]["candidate_no_upgrades_dir"]
        candidates = intersect_storage_candidates(data, no_upgrades_dir)

        # OPTIONAL OPTIONAL ADDITIONAL CAND STORAGE LOCATIONS
        if haskey(data["param"], "additional_cand")
            candidates = union!(candidates, string.(data["param"]["additional_cand"]))
        end

        all_busses = Set(keys(data["bus"]))
        non_candidates = setdiff(all_busses, union(candidates, nonzero_storage_nodes))
        non_candidates = Set(parse(Int, x) for x in non_candidates)

        for i in non_candidates
            fix(sigma[i], 0; force = true)
        end

        println("Number of candidates in Gurobi model: $(length(candidates))")
    end

    # Objective
    operational_weight = get(data["param"], "operational_weight", 1)
    @objective(model, Min,
        sum(s_power[i] * data["param"]["bess_power_cost"] + s_energy[i] * data["param"]["bess_energy_cost"] for i in 1:N) +
        sum(sigma[i] for i in 1:N) * data["param"]["storage_fixed_cost"] + 
        sum(data["param"]["cap_upgrade_cost"] * get_capacity_increment(data, a) * data["branch"]["$a"]["distance"] * gamma[a] for a in 1:E) +
        sum(
            data["param"]["representative_prob"][r] *
            (
                sum(
                    sum(compute_gen_cost(pg[r, g, t], data["gen"]["$g"]) for g in 1:G) +
                    sum(data["param"]["under_served_penalty"] * ue[r, i, t] for i in 1:N)
                for t in 1:T)
            )
        for r in 1:R) * operational_weight
    )
    
    return model
end

function compute_flows!(flows, pg_values, ue_values, ch_values, data, PTDF)
    # Preallocate net injection vector
    n_buses = length(data["bus"])
    net_injections = zeros(n_buses)  # PTDF matrix is (n_branches × n_buses)
    
    # For each timestep
    for r in 1:size(pg_values, 1), t in 1:size(pg_values, 3)
        fill!(net_injections, 0.0)

        # Compute net injections at each bus
        for i in 1:n_buses 
            # Sum generation at bus i
            net_injections[i] = sum(pg_values[r,g,t] 
                                    for g in 1:length(data["gen"]) 
                                    if data["gen"]["$g"]["gen_bus"] == i; 
                                    init=0.0)
            # Subtract load, add unserved energy, subtract charging
            net_injections[i] -= data["bus"]["$i"]["load"]["$r"][t]
            net_injections[i] += ue_values[r,i,t]
            net_injections[i] -= ch_values[r,i,t]
        end
        
        # Compute flows for all branches at once using PTDF
        flows[:,r,t] = PTDF * net_injections
    end
end

function save_tracked_constraints(simdir, model)
    # Convert dictionary to DataFrame
    df = DataFrame(arc = [k[1] for k in keys(model.ext[:tracked_constraints])],
    rep = [k[2] for k in keys(model.ext[:tracked_constraints])],
    time = [k[3] for k in keys(model.ext[:tracked_constraints])],
    tracked = collect(values(model.ext[:tracked_constraints])))

    # Save to CSV
    CSV.write(joinpath(simdir, "output", "tracked_constraints.csv"), df)
end

function save_power_injections(simdir, model, data)
    # Get solution values
    pg_values = value.(model[:pg])
    ue_values = value.(model[:ue])
    ch_values = value.(model[:ch])
    
    # Create storage for rows
    rows = []
    
    # Extract dimensions
    R = data["param"]["num_representatives"]
    G = length(data["gen"])
    N = length(data["bus"])
    T = data["param"]["num_hours"]
    
    # Gather all power injections
    for r in 1:R, t in 1:T
        # Generation
        for g in 1:G
            bus = data["gen"]["$g"]["gen_bus"]
            push!(rows, (rep=r, time=t, bus=bus, gen=g, variable="pg", value=pg_values[r, g, t]))
        end
        
        # Unserved energy
        for i in 1:N
            if ue_values[r, i, t] > 1e-6
                push!(rows, (rep=r, time=t, bus=i, gen=0, variable="ue", value=ue_values[r, i, t]))
            end
        end
        
        # Storage
        for i in 1:N
            if abs(ch_values[r, i, t]) > 1e-6
                push!(rows, (rep=r, time=t, bus=i, gen=0, variable="ch", value=ch_values[r, i, t]))
            end
        end
    end
    
    # Convert to DataFrame
    df = DataFrame(rows)
    
    # Save to CSV for each representative period
    for r in 1:R
        datestring = data["param"]["dates"][r]
        output_dir = joinpath(simdir, "output", datestring)
        mkpath(output_dir)
        
        df_rep = filter(row -> row.rep == r, df)
        CSV.write(joinpath(simdir, "output", "power_injections.csv"), df_rep)
    end
 end