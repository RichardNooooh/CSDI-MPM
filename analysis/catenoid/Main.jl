using ArgParse
include("../../src/ValidateArgs.jl")

function parseCommandLineArguments()
    s = ArgParseSettings(description="Runs a simulation of a catenoid.")
    @add_arg_table! s begin
        "--simulation_name"
            arg_type = String
            default = "Catenoid Membrane"
        "--mesh_file"
            help = "The file that defines the mesh for the simulation."
            arg_type = String
            default = "meshes/3D/catenoid_msh4.msh" # initially a cylinder
        "--grid_dims"
            help = "The lengths of the background grid along x, y, and z."
            nargs = 3
            arg_type = Float64
            default = [2.0, 2.0, 2.0]
        "--grid_num_cells"
            help = "Defines the number of cells along each axis."
            nargs = 3
            arg_type = Int64
            default = [8, 8, 8]
        "--dest"
            help = "The directory that holds the simulation results. This will hold the OUTPUT_NAME directory. " *
                   "If the directory does not exist at runtime, this program will crash."
            arg_type = String
            default = "output/analysis/catenoid/"
        "--output_name"
            help = "The name of the directory for the simulation results. The resulting path will " *
                   "then be DEST/OUTPUT_NAME/. If this directory does not exist at runtime, " *
                   "the program will crash. Also, note that paths in Julia depends on the directory " *
                   "the script is called from."
            arg_type = String
            default = "test/"
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
            default = 5000
        "--v_alpha"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC). Defaults to 1.0"
            arg_type = Float64
            default = 0.0
        "--area_density"
            help = "Defines the density of all surface material points in the simulation."
            arg_type = Float64
            default = 1.0
        "--surface_tension"
            help = "Defines the surface tension of the simulation."
            arg_type = Float64
            default = 1.0
        "--grid_r_min_active"
            help = "Defines the minimum radius of grid points that are 'active'. Grid points closer to the origin are ignored during computation."
            arg_type = Float64
            default = 0.99
        "--grid_r_max_active"
            help = "Defines the maximum radius of grid points that are 'active'. Grid points further from the origin are ignored during computation."
            arg_type = Float64
            default = 2.76
    end

    arguments = parse_args(s)
    validateMembraneArguments(arguments)
    return arguments
end

_arguments = parseCommandLineArguments()
include("Catenoid.jl")
run(_arguments)
