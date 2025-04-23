using ArgParse
include("../../src/ValidateArgs.jl")

function parseCommandLineArguments()
    s = ArgParseSettings(description="Runs a droplet simulation of a square droplet with a "*
                                     "incompressible viscous fluid domain. The domain has a symmetry boundary at x=0 and y=0.")
    @add_arg_table! s begin
        "--simulation_name"
            arg_type = String
            default = "Sessile Liquid Droplet"
        "--constitutive"
            arg_type = String
            default = "fluid"
        "--mesh_file"
            help = "The file that defines the mesh for the simulation."
            arg_type = String
            default = "meshes/2D/square_msh16.msh"
        "--grid_dims"
            help = "The lengths of the background grid along x and y."
            nargs = 2
            arg_type = Float64
            default = [1.0, 1.0]
        "--grid_num_cells"
            help = "Defines the number of cells along each axis."
            nargs = 2
            arg_type = Int64
            default = [32, 32]
        "--dest"
            help = "The directory that holds the simulation results. This will hold the OUTPUT_NAME directory. " *
                   "If the directory does not exist at runtime, this program will crash."
            arg_type = String
            default = "output/analysis_refactored/sessile_liquid_2d/"
        "--output_name"
            help = "The name of the directory for the simulation results. The resulting path will " *
                   "then be DEST/OUTPUT_NAME/. If this directory does not exist at runtime, " *
                   "the program will crash. Also, note that paths in Julia depends on the directory " *
                   "the script is called from."
            arg_type = String
            default = "test_h32_mp0.5/"
        "--t_init"
            help = "Defines the starting time of the simulation."
            arg_type = Float64
            default = 0.0
        "--t_delta"
            help = "Defines the time step in the explicit time integrations."
            arg_type = Float64
            default = 1.0e-3
        "--t_final"
            help = "Defines the ending time of the simulation. Simulation will not run longer than this."
            arg_type = Float64
            default = 1000.0
        "--t_load"
            help = "Defines the time in which gravity and surface tension will be " *
                   "at their correct magnitudes. Starting from t_init, the magnitudes of these " *
                   "loading forces will linearly increase from 0 to their unscaled magnitudes."
            arg_type = Float64
            default = 20.0
        "--t_record"
            help = "The simulation will record every \'t_record\' time steps."
            arg_type = Int64
            default = 1000
        "--v_alpha"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC). Defaults to 1.0"
            arg_type = Float64
            default = 0.0
        "--gravity", "-g"
            help = "Defines the magnitude of gravity in the -y direction."
            arg_type = Float64
            default = 0.0
        "--density"
            help = "Defines the density of all material points in the simulation."
            arg_type = Float64
            default = 1.0e-3 # g / mm^3
        "--bulk_modulus", "-K"
            help = "Defines the bulk modulus of the constitutive model."
            arg_type = Float64
            default = 0.0001 # MPa
        "--dynamic_viscosity"
            help = "Defines the dynamic viscosity of the constitutive model."
            arg_type = Float64
            default = 0.0 # MPa * ms
        "--surface_tension"
            help = "Defines the monolayer (droplet-air) surface tension of the simulation. This is " *
                   "only defined for the surface points on the curved surface."
            arg_type = Float64
            default = 72e-6 # N / mm
    end

    arguments = parse_args(s)
    validateArguments(arguments)
    return arguments
end

_arguments = parseCommandLineArguments()
include("SquareLiquidDroplet.jl")
run(_arguments)
