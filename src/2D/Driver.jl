include("../Util.jl")
include("ShapeFunctions.jl")
include("Grid.jl")
include("MaterialPoint.jl")
include("SurfacePoint.jl")
include("RigidPoint.jl")
include("ConstitutiveEquations.jl")
include("MPM.jl")
include("GmshVtkIO.jl")

"""
Initializes the grid with rollers on all 4 bounding planes of the grid,
adding extra cells on the border. `num_extra_cells` should be at least 1 to
prevent "Index out of bounds" errors at boundary.

# Arguments
- `arguments::Dict{String,Any}`: a dictionary of parameters
- `num_extra_cells::Int64`: number of extra grid cells on top and right boundaries

# Output
- Returns initialized `Grid` object
"""
function initializeGrid(arguments::Dict{String,Any}; num_extra_cells::Int64=1)::Grid
    @info("Initializing 2D $(arguments["simulation_name"]) Simulation")
    @info("Initializing background grid")

    grid_dims, grid_num_cells = arguments["grid_dims"], arguments["grid_num_cells"]
    h::Vector{Float64} = grid_dims ./ grid_num_cells
    grid::Grid = Grid(grid_dims[1] + num_extra_cells * h[1],
        grid_dims[2] + num_extra_cells * h[2],
        grid_num_cells[1] + 1 + num_extra_cells,
        grid_num_cells[2] + 1 + num_extra_cells)
    @info("Grid Setup:\n$(grid)")
    @threads for grid_point::GridPoint in grid.points
        # Zero-Friction Bounding Box
        if grid_point.position[1] < EPSILON || abs(grid_point.position[1] - grid.grid_length[1]) < EPSILON
            grid_point.is_fixed = Vec2{Bool}(true, grid_point.is_fixed[2])
        end
        if grid_point.position[2] < EPSILON || abs(grid_point.position[2] - grid.grid_length[2]) < EPSILON
            grid_point.is_fixed = Vec2{Bool}(grid_point.is_fixed[1], true)
        end
    end

    return grid
end


"""
Initializes the simulation with the material domain and surface domain.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `setBulkParameters!::Function`: a function that initializes the material domain with parameters. Takes in (material_domain, grid, arguments).
- `setSurfaceParameters!::Function`: a function that initializes the surface domain with parameters. Takes in (surface_domain, grid, arguments).
- `get_active_subset::Bool`: controls whether `active_surface_domain` is calculated

# Output
- Material domain
- Surface domain
- "Active" surface domain, the subset of the surface domain that will contribute a nonzero force to the grid nodes.
- 2-Tuple of VTKLookup tables for bulk and surface
"""
function initializeBodyMeshes(arguments::Dict{String,Any}, grid::Grid,
    setBulkParameters!::Function, setSurfaceParameters!::Function; get_active_subset::Bool=true)::Tuple{Vector{MaterialPoint},Vector{SurfacePoint},Vector{SurfacePoint},Tuple{VTKLookup,VTKLookup}}
    @info("Initializing bulk and surface from mesh file: $(arguments["mesh_file"])")
    material_domain::Vector{MaterialPoint}, surface_domain::Vector{SurfacePoint},
    bulk_lookup::VTKLookup, surf_lookup::VTKLookup = initialize2DGmshBulkAndSurface!(arguments["mesh_file"])
    setBulkParameters!(material_domain, grid, arguments)
    setSurfaceParameters!(surface_domain, grid, arguments)

    active_surface_domain::Vector{SurfacePoint} = surface_domain
    if get_active_subset
        active_surface_domain = Folds.collect(surface_point
                                              for surface_point::SurfacePoint in surface_domain
                                              if (surface_point.surface_tension > 0.0 || sum(surface_point.traction) > EPSILON))
    end

    return material_domain, surface_domain, active_surface_domain, (bulk_lookup, surf_lookup)
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
    rigid_vtk_lookup::VTKLookup = initialize2DGmshRigidBody!(arguments["rigid_mesh_file"])
    setRigidBodyParameters!(rigid_domain, grid, arguments)

    return rigid_domain, rigid_vtk_lookup
end


"""
Initializes the output folder directory, returning the directory names

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters

# Output
- Path names of newly created output subfolders
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
        mkdir(joinpath(vtk_output_folder, "bulk"))
        mkdir(joinpath(vtk_output_folder, "surface"))
    end

    analysis_output_folder = joinpath(destination_folder, "results/")
    if !isdir(analysis_output_folder)
        mkpath(analysis_output_folder)
    end

    return destination_folder, vtk_output_folder, analysis_output_folder
end


"""
Runs the main MPM loop with the initialized material/surface domains and grid.

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `grid::Grid`: the background grid
- `material_domain::Vector{MaterialPoint}`: a vector of `MaterialPoint`s
- `surface_domain::Vector{SurfacePoint}`: a vector of `SurfacePoint`s
- `active_surface_domain::Vector{SurfacePoint}`: a subset of `surface_domain` that will contribute nonzero forces to the grid
- `timeScaleFunctionBulk::Function`: a function that returns a value in [0, 1] to scale the gravity
- `timeScaleFunctionSurface::Function`: a function that returns a value in [0, 1] to scale the surface tension
- `recordVtkFiles!::Function`: a custom function that runs for every `t_record` time steps to record the VTK files
- `start_index::Int64` (optional keyword): the start index for the vtk files

# Output
- the time it took to run the algorithm in milliseconds
"""
function runMPM!(arguments::Dict{String,Any}, grid::Grid,
    material_domain::Vector{MaterialPoint}, surface_domain::Vector{SurfacePoint}, active_surface_domain::Vector{SurfacePoint},
    timeScaleFunctionBulk::Function, timeScaleFunctionSurface::Function,
    recordVtkFiles!::Function;
    start_index::Int64=Int64(arguments["t_init"] / arguments["t_delta"]))::Int64

    t_init::Float64, dt::Float64, t_final::Float64 = arguments["t_init"], arguments["t_delta"], arguments["t_final"]
    t_record::Int64 = arguments["t_record"]
    t_index::Int64 = start_index
    v_alpha::Float64 = arguments["v_alpha"]

    constitutiveEquation!::Function =
        arguments["constitutive"] == "hyperelastic" ? hyperelasticConstitutiveEquation! :
        arguments["constitutive"] == "fluid" ? fluidConstitutiveEquation! :
        error("Unknown constitutive model")

    start_time::DateTime = now()
    @showprogress 5 "SIMCODE | Running MPM Algorithm... " for t = t_init:dt:t_final
        # * plotting data
        if t_index % t_record == 0
            recordVtkFiles!(t, t_index)
        end
        t_index += 1

        t_scale_bulk::Float64 = timeScaleFunctionBulk(t)
        t_scale_surface::Float64 = timeScaleFunctionSurface(t)

        # Main MPM Logic
        resetGrid!(grid)
        materialToGrid!(grid, material_domain, t_scale_bulk)
        surfaceToGrid!(grid, active_surface_domain, t_scale_surface)
        updateGrid!(grid, dt)
        gridToMaterial!(material_domain, grid, dt, v_alpha, constitutiveEquation!)
        updateMaterialVertices!(material_domain, grid, dt)
        updateSurfaceVertices!(surface_domain, grid, dt)
    end
    end_time::DateTime = now()
    compute_time::Int64 = (end_time - start_time).value

    return compute_time
end


"""
Runs the main MPM loop with the initialized material/surface domains and grid.
Also applies dirichlet boundary conditions provided by rigid body domain

# Arguments
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `grid::Grid`: the background grid
- `material_domain::Vector{MaterialPoint}`: a vector of `MaterialPoint`s
- `surface_domain::Vector{SurfacePoint}`: a vector of `SurfacePoint`s
- `active_surface_domain::Vector{SurfacePoint}`: a subset of `surface_domain` that will contribute nonzero forces to the grid
- `rigid_domain::Vector{RigidPoint}`: a vector of `RigidPoint`s - same type of geometry as `MaterialPoint` (quads)
- `timeScaleFunctionBulk::Function`: a function that returns a value in [0, 1] to scale the gravity
- `timeScaleFunctionSurface::Function`: a function that returns a value in [0, 1] to scale the surface tension
- `recordVtkFiles!::Function`: a custom function that runs for every `t_record` time steps to record the VTK files
- `start_index::Int64` (optional keyword): the start index for the vtk files

# Output
- the time it took to run the algorithm in milliseconds
"""
function runMPM!(arguments::Dict{String,Any}, grid::Grid,
    material_domain::Vector{MaterialPoint}, surface_domain::Vector{SurfacePoint}, active_surface_domain::Vector{SurfacePoint},
    rigid_domain::Vector{RigidPoint},
    timeScaleFunctionBulk::Function, timeScaleFunctionSurface::Function,
    recordVtkFiles!::Function;
    start_index::Int64=Int64(arguments["t_init"] / arguments["t_delta"]))::Int64

    t_init::Float64, dt::Float64, t_final::Float64 = arguments["t_init"], arguments["t_delta"], arguments["t_final"]
    t_record::Int64 = arguments["t_record"]
    t_index::Int64 = start_index
    v_alpha::Float64 = arguments["v_alpha"]

    constitutiveEquation!::Function =
        arguments["constitutive"] == "hyperelastic" ? hyperelasticConstitutiveEquation! :
        arguments["constitutive"] == "fluid" ? fluidConstitutiveEquation! :
        error("Unknown constitutive model")

    start_time::DateTime = now()
    @showprogress 5 "SIMCODE | Running MPM Algorithm... " for t = t_init:dt:t_final
        # * plotting data
        if t_index % t_record == 0
            recordVtkFiles!(t, t_index)
        end
        t_index += 1

        t_scale_bulk::Float64 = timeScaleFunctionBulk(t)
        t_scale_surface::Float64 = timeScaleFunctionSurface(t)

        # Main MPM Logic
        resetGrid!(grid)
        materialToGrid!(grid, material_domain, t_scale_bulk)
        surfaceToGrid!(grid, active_surface_domain, t_scale_surface)
        updateRigidGrid!(grid, rigid_domain)    # rigid body DCs applies first, since `updateGrid!()` precalculates velocity
        updateGrid!(grid, dt)
        gridToMaterial!(material_domain, grid, dt, v_alpha, constitutiveEquation!)
        updateMaterialVertices!(material_domain, grid, dt)
        updateSurfaceVertices!(surface_domain, grid, dt)
        updateRigidVertices!(rigid_domain, grid, dt)
    end
    end_time::DateTime = now()
    compute_time::Int64 = (end_time - start_time).value

    return compute_time
end


"""
Records the setup for the material/surface domains, grid, and physical parameters to a file.

# Arguments
- `material_domain::Vector{MaterialPoint}`: a vector of `MaterialPoint`s
- `surface_domain::Vector{SurfacePoint}`: a vector of `SurfacePoint`s
- `active_surface_domain::Vector{SurfacePoint}`: a subset of `surface_domain` that will contribute nonzero forces to the grid
- `grid::Grid`: the background grid
- `arguments::Dict{String, Any}`: a dictionary of parameters
- `destination_folder::String`: directory to output setup file
- `extra_details::String`: string for additional details of setup

# Output
- Nothing
"""
function recordInitialSetupInfo(material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    active_surface_domain::Vector{SurfacePoint},
    grid::Grid, arguments::Dict{String,Any}, destination_folder::String, extra_details::String)::Nothing

    mass_volume::Vec2{Float64} = Folds.sum(Vec2{Float64}(material_point.mass, material_point.volume_init)
                                           for material_point::MaterialPoint in material_domain)
    mass::Float64, volume::Float64 = mass_volume
    total_area::Float64 = Folds.sum(surface_point.area for surface_point::SurfacePoint in surface_domain)
    total_active_area::Float64 = Folds.sum(surface_point.area for surface_point::SurfacePoint in active_surface_domain)

    setup_info = """
    $(arguments["simulation_name"]) Simulation

    Mesh File: \"$(arguments["mesh_file"])\"

    Initial configuration:
    Bulk:
        Total mass =   $mass
        Total volume = $volume
        Total number of material points = $(length(material_domain))
        Memory footprint: $(Base.summarysize(material_domain)) bytes
    Surface:
        Total area  =  $total_area
        Number of SP = $(length(surface_domain))
        Total "active" area = $total_active_area
        Number of active SP = $(length(active_surface_domain))
        Memory footprint: $(Base.summarysize(surface_domain)) bytes
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

    Gravity:           $(arguments["gravity"])
    Density:           $(arguments["density"])
    Velocity Alpha:    $(arguments["v_alpha"])

    Surface Tension:   $(arguments["surface_tension"])

    """

    if arguments["constitutive"] == "hyperelastic"
        E::Float64 = arguments["elastic_modulus"]
        v::Float64 = arguments["poisson_ratio"]
        λ::Float64, μ::Float64 = getLaméParameters(E, v)
        K::Float64 = E / (3 * (1 - 2 * v))
        setup_info *= """
        Hyperelastic Solid
            Elastic Modulus:   $(E)
            Poisson Ratio:     $(v)
            Bulk Modulus:      $(K)
            1st Lamé:          $(λ)
            2nd Lamé (Shear):  $(μ)
        """
    elseif arguments["constitutive"] == "fluid"
        setup_info *= """
        Viscous Fluid
            Bulk Modulus:      $(arguments["bulk_modulus"])
            Dynamic Viscosity: $(arguments["dynamic_viscosity"])
        """
    end

    setup_info *= "\n" * extra_details

    open(joinpath(destination_folder, ".parameters.txt"), "w") do file
        write(file, setup_info)
    end

    @info setup_info

    return nothing
end

