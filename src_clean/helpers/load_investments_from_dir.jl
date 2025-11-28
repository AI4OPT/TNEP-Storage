using CSV

function load_investments_from_dir(
    investment_dir::String,
    data::Dict;
    E::Int,
    N::Int,
    allow_missing::Bool=false
)
    """
    Load investment values from CSV files in a directory.
    
    Parameters:
    -----------
    investment_dir : String
        Directory containing investment CSV files
    data : Dict
        Model data dictionary (for storage scaling parameter)
    E : Int
        Number of lines (for validation and default zeros)
    N : Int
        Number of nodes (for validation and default zeros)
    allow_missing : Bool
        If true, return nothing for missing files; if false, return zeros
    
    Returns:
    --------
    (gamma, s_energy) : Tuple of (Union{Vector,Nothing}, Union{Vector,Nothing})
        Line and storage investment values
    """
    
    trans_file = joinpath(investment_dir, "line_investments.csv")
    storage_file = joinpath(investment_dir, "storage_investments.csv")
    
    # Initialize based on allow_missing flag
    gamma = allow_missing ? nothing : zeros(E)
    s_energy = allow_missing ? nothing : zeros(N)
    
    # Load line investments if file exists
    if isfile(trans_file)
        try
            trans_df = CSV.read(trans_file, DataFrame)
            gamma = trans_df[:, :Upgrade_Lvl]
            
            # Validate length
            if length(gamma) != E
                @warn "Line investments length ($(length(gamma))) doesn't match expected ($E). Using zeros."
                gamma = allow_missing ? nothing : zeros(E)
            else
                @info "Loaded $(length(gamma)) line investments from $trans_file"
            end
        catch e
            @warn "Failed to read line investments from $trans_file: $e"
            gamma = allow_missing ? nothing : zeros(E)
        end
    else
        @debug "Line investments file not found: $trans_file"
    end
    
    # Load storage investments if file exists
    if isfile(storage_file)
        try
            storage_df = CSV.read(storage_file, DataFrame)
            s_energy_raw = storage_df[:, :Storage_Energy]
            
            # Apply scaling if needed
            if get(data["param"], "storage_needs_scaling", false)
                s_energy = s_energy_raw .* data["param"]["storage_energy_size"]
            else
                s_energy = s_energy_raw
            end
            
            # Validate length
            if length(s_energy) != N
                @warn "Storage investments length ($(length(s_energy))) doesn't match expected ($N). Using zeros."
                s_energy = allow_missing ? nothing : zeros(N)
            else
                @info "Loaded $(length(s_energy)) storage investments from $storage_file"
            end
        catch e
            @warn "Failed to read storage investments from $storage_file: $e"
            s_energy = allow_missing ? nothing : zeros(N)
        end
    else
        @debug "Storage investments file not found: $storage_file"
    end
    
    return (gamma, s_energy)
end