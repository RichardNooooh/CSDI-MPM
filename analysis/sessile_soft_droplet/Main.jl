using ArgParse
include("../../src/ValidateArgs.jl")

function parseCommandLineArguments()
    s = ArgParseSettings(description="Runs a droplet simulation of a square droplet (No FBar) with a near-"*
                                     "incompressible viscous fluid domain. The domain has a symmetry boundary at x=0 and y=0.")
    @add_arg_table! s begin
        "--simulation_name"
            arg_type = String
            default = "Sessile Hyperelastic Droplet"
        "--constitutive"
            arg_type = String
            default = "hyperelastic"
        "--mesh_file"
            help = "The file that defines the mesh for the simulation."
            arg_type = String
            default = "meshes/2D/square_msh32.msh"
        "--grid_dims"
            help = "The lengths of the background grid along x and y."
            nargs = 2
            arg_type = Float64
            default = [1.0, 1.0]
        "--grid_num_cells"
            help = "Defines the number of cells along each axis."
            nargs = 2
            arg_type = Int64
            default = [64, 64]
        "--dest"
            help = "The directory that holds the simulation results. This will hold the OUTPUT_NAME directory. " *
                   "If the directory does not exist at runtime, this program will crash."
            arg_type = String
            default = "output/analysis/sessile_soft_2d/"
        "--output_name"
            help = "The name of the directory for the simulation results. The resulting path will " *
                   "then be DEST/OUTPUT_NAME/. If this directory does not exist at runtime, " *
                   "the program will crash. Also, note that paths in Julia depends on the directory " *
                   "the script is called from."
            arg_type = String
            default = "gamma_0.1/"
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
        "--elastic_modulus", "-E"
            help = "Defines the elastic_modulus of the constitutive model."
            arg_type = Float64
            default = 0.0001 # MPa
        "--poisson_ratio"
            help = "Defines the poisson_ratio of the constitutive model."
            arg_type = Float64
            default = 0.3 # MPa * ms
        "--surface_tension"
            help = "Defines the monolayer (droplet-air) surface tension of the simulation. This is " *
                   "only defined for the surface points on the curved surface."
            arg_type = Float64
            default = 72e-6 # N / mm
        "--surface_tension_factor"
            help = "Scales the surface tension parameter by this amount"
            arg_type = Float64
            default = 0.1
    end
    arguments = parse_args(s)
    validateArguments(arguments)
    return arguments
end

_arguments = parseCommandLineArguments()
include("SquareSoftDroplet.jl")
run(_arguments)
