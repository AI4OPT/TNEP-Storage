function recover_soc(simdir)
    data = JSON.parsefile(joinpath(simdir, "data.json"))
    filepath = joinpath(simdir, "output", "power_injections.csv")
    df = CSV.read(filepath, DataFrame)

    R = data["param"]["num_representatives"]
    T = data["param"]["num_hours"]
    η = data["param"]["bess_efficiency"]
    self_discharge = get(data["param"], "self_discharge", 0.0)
    storage_energy_size = data["param"]["storage_energy_size"]
    soc_init_ratio = get(data["param"], "soc_init_end_ratio", 0.5)

    # Load s_energy per bus from storage_investments.csv
    inv_path = joinpath(simdir, "output", "storage_investments.csv")
    inv_df = CSV.read(inv_path, DataFrame)
    s_energy = Dict(row.Node_Index => row.Storage_Energy for row in eachrow(inv_df))

    # Buses with installed storage
    storage_buses = sort([bus for (bus, e) in s_energy if e > 0])

    if isempty(storage_buses)
        @warn "No storage buses found in $inv_path"
        return nothing
    end

    # Build lookup: (rep, bus, time) => value
    ch_df = filter(r -> r.variable == "ch", df)
    dis_df = filter(r -> r.variable == "dis", df)
    ch_lookup = Dict((r.rep, r.bus, r.time) => r.value for r in eachrow(ch_df))
    dis_lookup = Dict((r.rep, r.bus, r.time) => r.value for r in eachrow(dis_df))

    # Build output: one row per (rep, bus), columns rep, bus, t_1..t_T
    n_rows = R * length(storage_buses)
    rep_col = Vector{Int}(undef, n_rows)
    bus_col = Vector{Int}(undef, n_rows)
    soc_mat = Matrix{Float64}(undef, n_rows, T)

    row_idx = 1
    for r in 1:R
        for bus in storage_buses
            e_cap = s_energy[bus] * storage_energy_size

            ch1 = get(ch_lookup, (r, bus, 1), 0.0)
            dis1 = get(dis_lookup, (r, bus, 1), 0.0)
            soc_mat[row_idx, 1] = soc_init_ratio * e_cap + ch1 * η - dis1 / η

            for t in 2:T
                ch_t = get(ch_lookup, (r, bus, t), 0.0)
                dis_t = get(dis_lookup, (r, bus, t), 0.0)
                soc_mat[row_idx, t] = (1 - self_discharge) * soc_mat[row_idx, t-1] + ch_t * η - dis_t / η
            end

            rep_col[row_idx] = r
            bus_col[row_idx] = bus
            row_idx += 1
        end
    end

    out_df = DataFrame(rep=rep_col, bus=bus_col)
    for t in 1:T
        out_df[!, Symbol("t_$t")] = soc_mat[:, t]
    end

    output_dir = joinpath(simdir, "output")
    mkpath(output_dir)
    out_path = joinpath(output_dir, "recovered_soc.csv")
    CSV.write(out_path, out_df)
    println("Wrote SOC recovery for $(length(storage_buses)) buses × $R reps to $out_path")

    return out_df
end
