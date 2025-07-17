using JSON
using Serialization

function save_master_problem(master, simdir; filename="master_problem.mps")
    """
    Save the master problem model to file for warm-starting
    """
    master_path = joinpath(simdir, "master_models")
    mkpath(master_path)
    
    # Save the MPS file
    mps_file = joinpath(master_path, filename)
    write_to_file(master, mps_file)
    
    # Save extension data (iterations tracking, etc.)
    extension_file = joinpath(master_path, replace(filename, ".mps" => "_extensions.jls"))
    extension_data = Dict(
        "y_raw" => master.ext[:y_raw],
        "y_eval" => master.ext[:y_eval],
        "y_core" => master.ext[:y_core],
        "iter" => master.ext[:iter],
        "gap" => master.ext[:gap],
        "stabilization_lambda" => master.ext[:stabilization_lambda],
        "stabilization_lambda_decr" => master.ext[:stabilization_lambda_decr],
        "INCIDENCE" => master.ext[:INCIDENCE]
    )
    
    serialize(extension_file, extension_data)
    
    println("[DEBUG] Master problem saved to $mps_file")
    println("[DEBUG] Extensions saved to $extension_file")
end

function load_master_problem(simdir, data; filename="master_problem.mps")
    """
    Load a previously saved master problem for warm-starting
    """
    master_path = joinpath(simdir, "master_models")
    mps_file = joinpath(master_path, filename)
    extension_file = joinpath(master_path, replace(filename, ".mps" => "_extensions.jls"))
    
    if !isfile(mps_file)
        println("[DEBUG] No saved master problem found at $mps_file")
        master, y, theta = define_master_ptdf(data)
        return master, y, theta
    end
    
    # Read the model from file (this creates a new model)
    optimizer = Gurobi.Optimizer
    master = read_from_file(mps_file)
    set_optimizer(master, optimizer)
    
    if !get(data["param"], "relaxed_first_stage", false)
        set_optimizer_attribute(master, "MIPGap", data["param"]["mip_gap"])
    end
    
    # Restore extension data
    if isfile(extension_file)
        extension_data = deserialize(extension_file)
        master.ext[:y_raw] = extension_data["y_raw"]
        master.ext[:y_eval] = extension_data["y_eval"]
        master.ext[:y_core] = extension_data["y_core"]
        master.ext[:iter] = extension_data["iter"]
        master.ext[:gap] = extension_data["gap"]
        master.ext[:stabilization_lambda] = extension_data["stabilization_lambda"]
        master.ext[:stabilization_lambda_decr] = extension_data["stabilization_lambda_decr"]
        master.ext[:INCIDENCE] = extension_data["INCIDENCE"]
        println("[DEBUG] Restored extension data with $(length(extension_data["y_raw"])) previous iterations")
    else
        # Initialize empty extensions if file doesn't exist
        master.ext[:y_raw] = Vector{Vector{Vector{Float64}}}()
        master.ext[:y_eval] = Vector{Vector{Vector{Float64}}}()
        master.ext[:y_core] = Vector{Vector{Vector{Float64}}}()
        master.ext[:gap] = Vector{Float64}()
        master.ext[:stabilization_lambda] = Vector{Float64}()
        master.ext[:stabilization_lambda_decr] = Vector{Float64}()
        master.ext[:INCIDENCE] = do_all_incidence(data)
    end
    
    # Reconstruct variable references based on the loaded model
    # This assumes the variable ordering is preserved
    E = length(data["branch"])
    N = length(data["bus"])
    
    # Get variable references from the loaded model
    all_vars = all_variables(master)
    
    # Reconstruct gamma variables (first E variables)
    gamma_vars = all_vars[1:E]
    master[:gamma] = gamma_vars
    
    # Reconstruct s_power variables (next N variables)
    s_power_vars = all_vars[E+1:E+N]
    master[:s_power] = s_power_vars
    
    # Reconstruct s_energy variables (next N variables)
    s_energy_vars = all_vars[E+N+1:E+2*N]
    master[:s_energy] = s_energy_vars
    
    # Reconstruct theta variable (last variable)
    theta_var = all_vars[end]
    master[:theta] = theta_var
    
    y = (master[:gamma], master[:s_power], master[:s_energy])
    
    println("[DEBUG] Successfully loaded master problem from $mps_file")
    return master, y, master[:theta]
end