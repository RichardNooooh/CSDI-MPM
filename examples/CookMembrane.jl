include("../src/ValidateArgs.jl")

_arguments = Dict{String,Any}([
    ("simulation_name", "Cook's Membrane"),
    ("constitutive", "hyperelastic"),
    ("mesh_file", "meshes/2D/cook_msh1.0.msh"),
    ("grid_dims", [80.0, 100.0]), # mm x mm
    ("grid_num_cells", [80, 100]),
    ("dest", "output/cook/"),
    ("output_name", "h1_ppc1/"),
    ("t_init", 0.0), # ms
    ("t_delta", 1e-3), # ms
    ("t_final", 1000.0), # ms
    ("t_load", 200),
    ("t_record", 1000),
    ("v_alpha", 0.0),
    ("elastic_modulus", 240.565),
    ("poisson_ratio", 0.3),
    ("surface_tension", 0.0),
    ("traction", 100.0 / 16.0),
    ("gravity", 0.0),
    ("density", 1.0e-3), # g / mm^3
])
validateArguments(_arguments)

include("../src/2D/Driver.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # *  Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments; num_extra_cells=1)
    for grid_point::GridPoint in grid.points
        if grid_point.position[1] < EPSILON || abs(grid_point.position[1] - grid.grid_length[1]) < EPSILON
            grid_point.is_fixed = Vec2{Bool}(true, true)
        end
    end

    material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    active_surface_domain::Vector{SurfacePoint},
    lookup_tables::Tuple{VTKLookup,VTKLookup} = initializeBodyMeshes(arguments, grid, setBulkParameters!, setSurfaceParameters!)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)

    recordInitialSetupInfo(material_domain, surface_domain, active_surface_domain,
        grid, arguments, destination_folder, extraSetupDetails(arguments))

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    bulk_data::BulkVTKData, pvd_bulk::WriteVTK.CollectionFile,
    vtk_bulk_cells::Vector{MeshCell} = initializeVTKBulkData(material_domain, destination_folder)
    surface_data::SurfaceVTKData, pvd_surf::WriteVTK.CollectionFile,
    vtk_surf_cells::Vector{MeshCell} = initializeVTKSurfaceData(surface_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)
    writeVTIGrid(grid, destination_folder)

    t_record::Int64 = arguments["t_record"]
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKBulk(material_domain, bulk_data, lookup_tables[1], vtk_bulk_cells, pvd_bulk,
            t_index, t_record, t, joinpath(vtk_output_folder, "bulk"))
        writeVTKSurface(surface_domain, surface_data, lookup_tables[2], vtk_surf_cells, pvd_surf,
            t_index, t_record, t, joinpath(vtk_output_folder, "surface"))
        writeVTKGrid(grid, grid_data, pvd_grid,
            t_index, t_record, t, joinpath(vtk_output_folder, "grid"))
    end

    # * ======================================
    # * Time Scale Function
    # * ======================================
    t_load::Float64 = arguments["t_load"]
    function timeScaleFunction(t::Float64)::Float64
        return minimum((t_load, t)) / t_load
    end

    # * ======================================
    # * Running MPM Algorithm
    # * ======================================
    compute_time::Int64 = runMPM!(arguments, grid, material_domain, surface_domain, active_surface_domain,
        timeScaleFunction, timeScaleFunction, recordVtkFiles!)

    # * ======================================
    # * Final Output Files
    # * ======================================
    @info("Saving PVD files")
    vtk_save(pvd_bulk)
    vtk_save(pvd_surf)
    vtk_save(pvd_grid)

    @info("Recording compute time")
    open(analysis_output_folder * "compute_time_ms.txt", "w") do file
        write(file, "$compute_time\n")
    end

    @info("Recording maximum grid velocity")
    max_grid_speed::Float64 = Folds.maximum(norm(grid_point.velocity_next)
                                            for grid_point::GridPoint in grid.points)
    open(joinpath(analysis_output_folder, "max_grid_speed.txt"), "w") do file
        write(file, "$max_grid_speed\n")
    end

    @info("Recording All (Active) Surface Point Positions")
    open(joinpath(analysis_output_folder, "vertex_list.tsv"), "w") do file
        write(file, "x1\ty1\tx2\ty2\n")
        for surface_point::SurfacePoint in active_surface_domain
            p1::Vec2{Float64}, p2::Vec2{Float64} = surface_point.vertices
            write(file, "$(p1[1])\t$(p1[2])\t")
            write(file, "$(p2[1])\t$(p2[2])\n")
        end
    end

    @info("Recording All Material Point Pressures and Volumes")
    open(joinpath(analysis_output_folder, "pressure_list.tsv"), "w") do file
        write(file, "pressure\tvolume\n")
        for material_point::MaterialPoint in material_domain
            write(file, "$(material_point.pressure)\t$(material_point.volume)\n")
        end
    end

    @info("Recording total simulation time")
    open(joinpath(analysis_output_folder, "compute_time_ms.txt"), "w") do file
        write(file, "$compute_time\n")
    end

    @info("Finished simulation! :)")
    open(joinpath(analysis_output_folder, "COMPLETE"), "w") do file
        write(file, "This is the receipt that the simulation is finished.\n")
    end
    return nothing
end


function setBulkParameters!(material_domain::Vector{MaterialPoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    ρ::Float64 = arguments["density"]
    E::Float64 = arguments["elastic_modulus"]
    v::Float64 = arguments["poisson_ratio"]

    λ::Float64, μ::Float64 = getLaméParameters(E, v)

    @threads for material_point::MaterialPoint in material_domain
        m::Float64 = material_point.volume_init * ρ

        material_point.mass = m
        material_point.density = ρ
        material_point.parameter_1 = λ
        material_point.parameter_2 = μ
        material_point.force_external = ZERO_VEC2
        material_point.velocity = ZERO_VEC2

        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices,
            material_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function setSurfaceParameters!(surface_domain::Vector{SurfacePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    T::Float64 = arguments["traction"]

    for surface_point::SurfacePoint in surface_domain
        surface_point.surface_tension = 0.0
        surface_point.traction = ZERO_VEC2
        # Right-hand Traction Boundary
        if abs((getSumOfVertexInDimension(surface_point, 1) / 2) - 48.0) < EPSILON
            surface_point.traction = Vec2{Float64}(0.0, T)
        else
            surface_point.traction = ZERO_VEC2
        end

        surface_point.connected_grid_array_length = getAllConnectedGrid!(surface_point.connected_grid_indices,
            surface_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function extraSetupDetails(arguments::Dict{String,Any})::String
    return """
    Notes
    - Units: grams (g), millimeters (mm), milliseconds (ms). Pressure is in MPa, force is in N.
    - Traction: $(arguments["traction"])    
    - Traction is loaded linearly over $(arguments["t_load"])
    """
end


run(_arguments)
