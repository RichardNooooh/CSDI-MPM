using ArgParse
include("../../src/ValidateArgs.jl")

function parseCommandLineArguments()
    s = ArgParseSettings(description="Runs a droplet simulation of a 2 droplets pulled apart")
    @add_arg_table! s begin
        "--simulation_name"
            arg_type = String
            default = "Droplet Pull"
        "--constitutive"
            arg_type = String
            default = "hyperelastic"
        "--mesh_file"
            help = "The file that defines the mesh for the simulation."
            arg_type = String
            default = "meshes/2D/twodisk_r1.0_r1.0_xa1.25_xb4.25_msh16.msh"
        "--rigid_mesh_file"
            help = "The file that defines the mesh for the simulation."
            arg_type = String
            default = "meshes/2D/twodisk_r0.25_r0.25_xa1.25_xb4.25_msh16.msh"
        "--grid_dims"
            help = "The lengths of the background grid along x and y."
            nargs = 2
            arg_type = Float64
            default = [5.5, 3.0]
        "--grid_num_cells"
            help = "Defines the number of cells along each axis."
            nargs = 2
            arg_type = Int64
            default = [88, 48]
        "--dest"
            help = "The directory that holds the simulation results. This will hold the OUTPUT_NAME directory. " *
                   "If the directory does not exist at runtime, this program will crash."
            arg_type = String
            default = "output/analysis/droplet_pull/"
        "--output_name"
            help = "The name of the directory for the simulation results. The resulting path will " *
                   "then be DEST/OUTPUT_NAME/. If this directory does not exist at runtime, " *
                   "the program will crash. Also, note that paths in Julia depends on the directory " *
                   "the script is called from."
            arg_type = String
            default = "no_st/"
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
            default = 500.0
        "--t_load"
            help = "Defines the time in which gravity and surface tension will be " *
                   "at their correct magnitudes. Starting from t_init, the magnitudes of these " *
                   "loading forces will linearly increase from 0 to their unscaled magnitudes."
            arg_type = Float64
            default = 20.0
        "--t_record"
            help = "The simulation will record VTK files every \'t_record\' time steps."
            arg_type = Int64
            default = 1000
        "--v_alpha"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC) for the second pull-apart stage"
            arg_type = Float64
            default = 1.0
        "--gravity", "-g"
            help = "Defines the magnitude of gravity in the -y direction."
            arg_type = Float64
            default = 0.0
        "--density"
            help = "Defines the density of all material points in the simulation."
            arg_type = Float64
            default = 1.0e-3
        "--elastic_modulus", "-E"
            help = "Defines the elastic modulus of the hyperelastic block."
            arg_type = Float64
            default = 0.0001 # N / mm^2 == MPa
        "--poisson_ratio"
            help = "Defines the dynamic viscosity of the hyperelastic block."
            arg_type = Float64
            default = 0.3
        "--surface_tension"
            help = "Defines the monolayer (droplet-air) surface tension of the simulation."
            arg_type = Float64
            default = 0.0
        "--t_velocity_start"
            help = "Defines start time of rigid body pull"
            arg_type = Float64
            default = 200.0
        "--velocity"
            help = "Defines how fast the rigid body will be moving away from each other."
            arg_type = Float64
            default = 0.001
        "--v_alpha_initial"
            help = "Defines velocity update (alpha = 1.0 for FLIP, 0.0 for PIC) for the first compressive stage"
            arg_type = Float64
            default = 0.0
        "--translation_vector_1"
            help = "Shifts left droplet by this vector. Center at (1.25, 1.25)."
            nargs = 2
            arg_type = Float64
            default = [0.5-1/128, 0.25]
        "--translation_vector_2"
            help = "Shifts right droplet by this vector. Center at (4.25, 1.25)."
            nargs = 2
            arg_type = Float64
            default = [-0.5+1/128, 0.25]
    end

    arguments = parse_args(s)
    validateArguments(arguments)
    return arguments
end

_arguments = parseCommandLineArguments()
include("DropletPullStudy.jl")
run(_arguments)
