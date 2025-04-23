using ArgParse
include("../../src/ValidateArgs.jl")

function parseCommandLineArguments()
    s = ArgParseSettings(description="Runs a pendant droplet simulation.")
    @add_arg_table! s begin
        "--simulation_name"
            arg_type = String
            default = "Pendant Droplet"
        "--constitutive"
            arg_type = String
            default = "fluid"
        "--mesh_file"
            help = "The file that defines the hemisphere mesh for the simulation."
            arg_type = String
            default = "meshes/3D/pendant_msh32.msh"
        "--grid_dims"
            help = "The lengths of the background grid along x, y, and z."
            nargs = 3
            arg_type = Float64
            default = [0.75, 0.75, 3.0]
        "--grid_num_cells"
            help = "Defines the number of cells along each axis."
            nargs = 3
            arg_type = Int64
            default = [24, 24, 96]
        "--translation_vector"
            help = "Shifts mesh by this vector."
            nargs = 3
            arg_type = Float64
            default = [0.0, 0.0, 0.0]
        "--dest"
            help = "The directory that holds the simulation results. This will hold the OUTPUT_NAME directory. " *
                   "If the directory does not exist at runtime, this program will crash."
            arg_type = String
            default = "output/analysis/pendant/"
        "--output_name"
            help = "The name of the directory for the simulation results. The resulting path will " *
                   "then be DEST/OUTPUT_NAME/. If this directory does not exist at runtime, " *
                   "the program will crash. Also, note that paths in Julia depends on the directory " *
                   "the script is called from."
            arg_type = String
            default = "h32_msh32/"
        "--t_init"
            help = "Defines the starting time of the simulation."
            arg_type = Float64
            default = 0.0
        "--t_delta"
            help = "Defines the time step in the explicit time integrations."
            arg_type = Float64
            default = 1e-3
        "--t_final"
            help = "Defines the ending time of the simulation."
            arg_type = Float64
            default = 600.0
        "--t_record"
            help = "The simulation will record every \'t_record\' time steps."
            arg_type = Int64
            default = 1000
        "--t_load_g_start"
            help = "Time to start gravity"
            arg_type = Float64
            default = 50.0
        "--t_load_g_max"
            help = "Time to gravity max from t_load_g_start"
            arg_type = Float64
            default = 50.0
        "--t_load_st_max"
            help = "Time to surface tension max"
            arg_type = Float64
            default = 50.0
        "--t_load_st_decay_start"
            help = "Time to surface tension decay start"
            arg_type = Float64
            default = 150.0
        "--t_load_st_decay_factor"
            help = "Time to surface tension decay factor"
            arg_type = Float64
            default = 0.005
        "--v_alpha"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC) for 2nd stage"
            arg_type = Float64
            default = 1.0
        "--v_alpha_initial"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC) for 1st stage"
            arg_type = Float64
            default = 0.0
        "--gravity", "-g"
            help = "Defines the magnitude of gravity in the -z direction."
            arg_type = Float64
            default = 9.8e-3
        "--density"
            help = "Defines the density of all material points in the simulation."
            arg_type = Float64
            default = 1e-3
        "--bulk_modulus", "-K"
            help = "Defines the bulk modulus of the constitutive model."
            arg_type = Float64
            default = 0.0004 # MPa
        "--dynamic_viscosity"
            help = "Defines the dynamic viscosity of the constitutive model."
            arg_type = Float64
            default = 1e-6 # MPa * ms
        "--surface_tension"
            help = "Defines the monolayer (droplet-air) surface tension of the simulation"
            arg_type = Float64
            default = 50e-6
    end

    arguments = parse_args(s)
    validateArguments(arguments)
    return arguments
end

_arguments = parseCommandLineArguments()
include("PendantDroplet.jl")
run(_arguments)
