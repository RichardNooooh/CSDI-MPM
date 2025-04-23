include("../Util.jl")
include("ShapeFunctions.jl")
include("Grid.jl")
include("MembranePoint.jl")
include("RigidPoint.jl")
include("MPM.jl")
include("GmshVtkIO.jl")

"""
Initializes the grid with rollers on all 6 bounding planes of the grid,
adding extra cells on the border
"""
function initializeGrid(arguments::Dict{String,Any}; num_extra_cells::Int64=1)::Grid
    @info("Initializing 3D $(arguments["simulation_name"]) Simulation")
    @info("Initializing background grid")

    grid_dims, grid_num_cells = arguments["grid_dims"], arguments["grid_num_cells"]
    h::Vector{Float64} = grid_dims ./ grid_num_cells
    grid::Grid = Grid(grid_dims[1] + num_extra_cells * h[1],
        grid_dims[2] + num_extra_cells * h[2],
        grid_dims[3] + num_extra_cells * h[3],
        grid_num_cells[1] + 1 + num_extra_cells,
        grid_num_cells[2] + 1 + num_extra_cells,
        grid_num_cells[3] + 1 + num_extra_cells)
    @info("Grid Setup:\n$(grid)")
    @threads for grid_point::GridPoint in grid.points
        # Symmetry boundary conditions
        if grid_point.position[1] < EPSILON || abs(grid_point.position[1] - grid.grid_length[1]) < EPSILON
            grid_point.is_fixed = Vec3{Bool}(true, grid_point.is_fixed[2], grid_point.is_fixed[3])
        end
        if grid_point.position[2] < EPSILON || abs(grid_point.position[2] - grid.grid_length[2]) < EPSILON
            grid_point.is_fixed = Vec3{Bool}(grid_point.is_fixed[1], true, grid_point.is_fixed[3])
        end
        if grid_point.position[3] < EPSILON || abs(grid_point.position[3] - grid.grid_length[3]) < EPSILON
            grid_point.is_fixed = Vec3{Bool}(grid_point.is_fixed[1], grid_point.is_fixed[2], true)
        end
    end

    return grid
end


"""
Initializes the simulation with the membrane domain.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `setMembraneParameters!::Function`: a function that initializes the membrane domain with parameters. Takes in (membrane_domain, grid, arguments).

# Output
- Membrane domain
- VTKLookup table
"""
function initializeBodyMesh(arguments::Dict{String,Any}, grid::Grid,
    setMembraneParameters!::Function)::Tuple{Vector{MembranePoint},VTKLookup}
    @info("Initializing membrane from mesh file: $(arguments["mesh_file"])")
    membrane_domain::Vector{MembranePoint},
    membrane_lookup::VTKLookup = initialize3DGmshSurface!(arguments["mesh_file"])
    setMembraneParameters!(membrane_domain, grid, arguments)

    return membrane_domain, membrane_lookup
end

"""
Initializes the simulation with the rigid domain.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `setRigidBodyParameters!::Function`: a function that initializes the rigid body domain with parameters. Takes in (rigid_domain, grid).

# Output
- Rigid body domain
- VTKLookup table
"""
function initializeRigidMesh(arguments::Dict{String,Any}, grid::Grid,
    setRigidBodyParameters!::Function)::Tuple{Vector{RigidPoint},VTKLookup}
    @info("Initializing rigid body from mesh file: $(arguments["rigid_mesh_file"])")
    rigid_domain::Vector{RigidPoint},
    rigid_vtk_lookup::VTKLookup = initialize3DGmshRigidSurface!(arguments["rigid_mesh_file"])
    setRigidBodyParameters!(rigid_domain, grid)

    return rigid_domain, rigid_vtk_lookup
end


"""
Initializes the output folder directory, returning the directory names
"""
function initializeOutputDirectory(arguments::Dict{String,Any})::Tuple{String,String,String}
    @info("Creating output file structure in $(arguments["dest"])")
    destination_folder::String = joinpath(arguments["dest"], arguments["output_name"])
    if !isdir(destination_folder)
        mkpath(destination_folder)
    end

    vtk_output_folder = joinpath(destination_folder, "vtk/")
    if !isdir(vtk_output_folder)
        mkpath(vtk_output_folder)
        mkdir(joinpath(vtk_output_folder, "grid"))
        mkdir(joinpath(vtk_output_folder, "membrane"))
    end

    analysis_output_folder = joinpath(destination_folder, "results/")
    if !isdir(analysis_output_folder)
        mkpath(analysis_output_folder)
    end

    return destination_folder, vtk_output_folder, analysis_output_folder
end


"""
Runs the main MPM loop with the initialized membrane domain and grid.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `grid::Grid`: the background grid
- `membrane_domain::Vector{MembranePoint}`: a vector of `MembranePoint`s
- `timeScaleFunction::Function`: a function that returns a value in [0, 1] to scale the surface tension
- `recordVtkFiles!::Function`: a custom function that runs for every `t_record` time steps to record the VTK files
- `start_index::Int64` (optional keyword): the start index for the vtk files

# Output
- the time it took to run the algorithm in milliseconds
"""
function runMPM!(arguments::Dict{String,Any}, grid::Grid,
    membrane_domain::Vector{MembranePoint},
    timeScaleFunction::Function,
    recordVtkFiles!::Function;
    start_index::Int64=Int64(arguments["t_init"] / arguments["t_delta"]))::Int64
    t_init::Float64, dt::Float64, t_final::Float64 = arguments["t_init"], arguments["t_delta"], arguments["t_final"]
    t_record::Int64 = arguments["t_record"]

    t_index::Int64 = start_index
    v_alpha::Float64 = arguments["v_alpha"]

    start_time::DateTime = now()
    @showprogress 5 "SIMCODE | Running MPM Algorithm... " for t = t_init:dt:t_final
        if t_index % t_record == 0
            recordVtkFiles!(t, t_index)
        end
        t_index += 1
        t_scale::Float64 = timeScaleFunction(t)

        # Main MPM Logic
        resetGrid!(grid.points)
        membraneToGrid!(grid, membrane_domain, t_scale)
        updateGrid!(grid.points, dt)
        gridToMembrane!(membrane_domain, grid, v_alpha)
        updateMembraneVertices!(membrane_domain, grid, dt)
    end
    end_time::DateTime = now()
    compute_time::Int64 = (end_time - start_time).value

    return compute_time
end

"""
Ensures that all grid points are zero'd-out
"""
function gridSanityCheck(grid_points::Vector{GridPoint})::Nothing
    badActiveGrid::Threads.Atomic{Bool} = Threads.Atomic{Bool}(false)
    @threads for grid_point::GridPoint in grid_points
        if badActiveGrid[]
            continue
        end
        if grid_point.mass > EPSILON || sum(grid_point.momentum) > EPSILON || sum(grid_point.force) > EPSILON
            badActiveGrid[] = true
        end
    end

    if badActiveGrid[]
        println("SIMCODE    | !!!! A grid point with mass > EPSILON was found. Your active grid is not large enough or the mesh exploded")
        for grid_point::GridPoint in grid_points
            if grid_point.mass > EPSILON || sum(grid_point.momentum) > EPSILON || sum(grid_point.force) > EPSILON
                println("SIMCODE    | Bad grid node at $(grid_point.position) with mass = $(grid_point.mass), momentum = $(grid_point.momentum), force = $(grid_point.force).")
                break
            end
        end
        exit(1)
    end
end


"""
Runs the main MPM loop with the initialized membrane domain and grid. Utilizes a subset of the grid points for resetting and updating.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `grid::Grid`: the background grid
- `active_grid_points::Vector{GridPoint}`: a subset of grid points that are used by the simulation. A sanity check function is called every `t_record` time steps just in case.
- `membrane_domain::Vector{MembranePoint}`: a vector of `MembranePoint`s
- `timeScaleFunction::Function`: a function that returns a value in [0, 1] to scale the surface tension
- `recordVtkFiles!::Function`: a custom function that runs for every `t_record` time steps to record the VTK files
- `start_index::Int64` (optional keyword): the start index for the vtk files

# Output
- the time it took to run the algorithm in milliseconds
"""
function runMPM!(arguments::Dict{String,Any}, grid::Grid, active_grid_points::Vector{GridPoint},
    membrane_domain::Vector{MembranePoint},
    timeScaleFunction::Function,
    recordVtkFiles!::Function;
    start_index::Int64=Int64(arguments["t_init"] / arguments["t_delta"]))::Int64
    @info("Shrunk active grid from $(length(grid.points)) grid points to $(length(active_grid_points))")

    t_init::Float64, dt::Float64, t_final::Float64 = arguments["t_init"], arguments["t_delta"], arguments["t_final"]
    t_record::Int64 = arguments["t_record"]

    t_index::Int64 = start_index
    v_alpha::Float64 = arguments["v_alpha"]

    start_time::DateTime = now()
    @showprogress 5 "SIMCODE | Running MPM Algorithm... " for t = t_init:dt:t_final
        if t_index % t_record == 0
            recordVtkFiles!(t, t_index)
        end
        t_index += 1
        t_scale::Float64 = timeScaleFunction(t)

        # Main MPM Logic
        resetGrid!(active_grid_points)
        if t_index % t_record == 0
            gridSanityCheck(grid.points)
        end

        membraneToGrid!(grid, membrane_domain, t_scale)
        updateGrid!(active_grid_points, dt)
        gridToMembrane!(membrane_domain, grid, v_alpha)
        updateMembraneVertices!(membrane_domain, grid, dt)

    end
    end_time::DateTime = now()
    compute_time::Int64 = (end_time - start_time).value

    return compute_time
end


"""
Runs the main MPM loop with the initialized membrane domain and grid. Utilizes a subset of the grid points for resetting and updating.

Also uses fixed rigid points for boundary conditions.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `grid::Grid`: the background grid
- `active_grid_points::Vector{GridPoint}`: a subset of grid points that are used by the simulation. A sanity check function is called every `t_record` time steps just in case.
- `active_fixed_grid_points::Vector{GridPoint}`: a subset of `active_grid_points` that are in contact with a fixed rigid body
- `membrane_domain::Vector{MembranePoint}`: a vector of `MembranePoint`s
- `timeScaleFunction::Function`: a function that returns a value in [0, 1] to scale the surface tension
- `recordVtkFiles!::Function`: a custom function that runs for every `t_record` time steps to record the VTK files
- `start_index::Int64` (optional keyword): the start index for the vtk files

# Output
- the time it took to run the algorithm in milliseconds
"""
function runMPM!(arguments::Dict{String,Any}, grid::Grid,
    active_grid_points::Vector{GridPoint}, active_fixed_rigid_grid::Vector{GridPoint},
    membrane_domain::Vector{MembranePoint},
    timeScaleFunction::Function,
    recordVtkFiles!::Function;
    start_index::Int64=Int64(arguments["t_init"] / arguments["t_delta"]))::Int64
    @info("Shrunk active grid from $(length(grid.points)) grid points to $(length(active_grid_points))")
    @info("Fixed active grid from $(length(active_grid_points)) to $(length(active_fixed_rigid_grid))")

    t_init::Float64, dt::Float64, t_final::Float64 = arguments["t_init"], arguments["t_delta"], arguments["t_final"]
    t_record::Int64 = arguments["t_record"]

    t_index::Int64 = start_index
    v_alpha::Float64 = arguments["v_alpha"]

    start_time::DateTime = now()
    @showprogress 5 "SIMCODE | Running MPM Algorithm... " for t = t_init:dt:t_final
        if t_index % t_record == 0
            recordVtkFiles!(t, t_index)
        end
        t_index += 1
        t_scale::Float64 = timeScaleFunction(t)

        # Main MPM Logic
        resetGrid!(active_grid_points)
        if t_index % t_record == 0
            gridSanityCheck(grid.points)
        end

        membraneToGrid!(grid, membrane_domain, t_scale)

        updateRigidGridFixed!(active_fixed_rigid_grid)
        updateGrid!(active_grid_points, dt)

        gridToMembrane!(membrane_domain, grid, v_alpha)
        updateMembraneVertices!(membrane_domain, grid, dt)
    end
    end_time::DateTime = now()
    compute_time::Int64 = (end_time - start_time).value

    return compute_time
end


"""
Records the setup for the membrane domain, grid, and physical parameters to a file.
"""
function recordInitialSetupInfo(membrane_domain::Vector{MembranePoint},
    grid::Grid, arguments::Dict{String,Any}, destination_folder::String, extra_details::String)::Nothing
    mass_area::Vec2{Float64} = Folds.sum(Vec2{Float64}(membrane_point.mass, membrane_point.area)
                                         for membrane_point::MembranePoint in membrane_domain)
    mass::Float64, area::Float64 = mass_area

    setup_info = """
    $(arguments["simulation_name"]) Simulation

    Mesh File: \"$(arguments["mesh_file"])\"

    Initial configuration:
    Membrane:
        Total mass = $mass
        Total area  =  $area
        Number of MP = $(length(membrane_domain))
        Memory footprint: $(Base.summarysize(membrane_domain)) bytes
    Grid:
        Grid length =      $(grid.grid_length)
        Grid node count =  $(grid.num_nodes)
        Grid cell length = $(grid.cell_length)
        Memory footprint: $(Base.summarysize(grid)) bytes

    Threads: $(Threads.nthreads())

    Initial Time:      $(arguments["t_init"])
    Time Step:         $(arguments["t_delta"])
    Final Time:        $(arguments["t_final"])
    Records Every      $(arguments["t_record"]) Time Steps
        For a total of $(Int64(ceil((arguments["t_final"] - arguments["t_init"]) / arguments["t_delta"] / arguments["t_record"]))) Recorded Time Steps

    Density:           $(arguments["area_density"])
    Velocity Alpha:    $(arguments["v_alpha"])

    Surface Tension:   $(arguments["surface_tension"])
    """

    setup_info *= "\n" * extra_details

    open(joinpath(destination_folder, ".parameters.txt"), "w") do file
        write(file, setup_info)
    end

    @info setup_info

    return nothing
end

