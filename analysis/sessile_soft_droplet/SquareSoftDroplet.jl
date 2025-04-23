include("../../src/2D/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # *  Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments; num_extra_cells=1)

    material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    active_surface_domain::Vector{SurfacePoint},
    lookup_tables::Tuple{VTKLookup,VTKLookup} = initializeBodyMeshes(arguments, grid, setBulkParameters!, setSurfaceParameters!)

    tracked_material_point::MaterialPoint, tracked_vertex::Int64 = trackFarthestCorner(material_domain)
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

    @info("Record elastocapillary number of the simulation")
    open(joinpath(analysis_output_folder, "elastocapnum.txt"), "w") do file
        write(file, "$(getElastocapillaryNumber(tracked_material_point, tracked_vertex, arguments))\n")
    end

    @info("Finished simulation! :)")
    open(joinpath(analysis_output_folder, "COMPLETE"), "w") do file
        write(file, "This is the receipt that the simulation is finished.\n")
    end
    return nothing
end


function setBulkParameters!(material_domain::Vector{MaterialPoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    E::Float64 = arguments["elastic_modulus"]
    v::Float64 = arguments["poisson_ratio"]
    g::Float64 = arguments["gravity"]
    ρ::Float64 = arguments["density"]

    λ::Float64, μ::Float64 = getLaméParameters(E, v)

    @threads for material_point::MaterialPoint in material_domain
        m::Float64 = material_point.volume_init * ρ

        material_point.mass = m
        material_point.density = ρ
        material_point.parameter_1 = λ
        material_point.parameter_2 = μ
        material_point.force_external = Vec2{Float64}(0.0, -m * g)
        material_point.velocity = ZERO_VEC2

        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices,
            material_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function setSurfaceParameters!(surface_domain::Vector{SurfacePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    γ::Float64 = arguments["surface_tension"] * arguments["surface_tension_factor"]

    @threads for surface_point::SurfacePoint in surface_domain
        # Symmetry Lines x = 0 and y = 0
        if getSumOfVertexInDimension(surface_point, 1) < EPSILON || getSumOfVertexInDimension(surface_point, 2) < EPSILON
            surface_point.surface_tension = 0.0
            # Curved Droplet-Air Surface
        else
            surface_point.surface_tension = γ
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
    - Surface tension is loaded linearly over $(arguments["t_load"])
    - Above surface tension is scaled by a constant factor of $(arguments["surface_tension_factor"])
    - Symmetry boundaries on x=0 and y=0
    """
end
