"""
Checks that all required keys exist in `arguments` and that the key 
"constitutive" has a valid value ("hyperelastic" or "fluid"). Throws an 
error if any requirement is not met.
"""
function validateArguments(arguments::Dict{String, Any})::Nothing
    required_keys = [
        "simulation_name",
        "grid_dims",
        "grid_num_cells",
        "mesh_file",
        "dest",
        "output_name",
        "t_init",
        "t_delta",
        "t_final",
        "t_record",
        "v_alpha",
        "constitutive"
    ]

    for key in required_keys
        if !haskey(arguments, key)
            error("Missing required argument: $key")
        end
    end

    valid_constitutive = ["hyperelastic", "fluid"]
    if !(arguments["constitutive"] in valid_constitutive)
        error("Invalid value for key 'constitutive': $(arguments["constitutive"]). Must be one of: $(valid_constitutive).")
    end

    return nothing
end

"""
Checks that all required keys exist in `arguments` and that the key 
"constitutive" has a valid value ("hyperelastic" or "fluid"). Throws an 
error if any requirement is not met.
"""
function validateMembraneArguments(arguments::Dict{String, Any})::Nothing
    required_keys = [
        "simulation_name",
        "grid_dims",
        "grid_num_cells",
        "mesh_file",
        "dest",
        "output_name",
        "t_init",
        "t_delta",
        "t_final",
        "t_record",
        "v_alpha"
    ]

    for key in required_keys
        if !haskey(arguments, key)
            error("Missing required argument: $key")
        end
    end

    return nothing
end

