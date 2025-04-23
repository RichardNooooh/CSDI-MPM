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

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)

    expected_radius::Float64, expected_pressure::Float64 = calculateExpectedRadiiAndPressure(arguments)

    recordInitialSetupInfo(material_domain, surface_domain, active_surface_domain,
        grid, arguments, destination_folder, extraSetupDetails(arguments, expected_radius, expected_pressure))

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

    numMP::Int64 = length(material_domain)
    pressure_diff::Array{Float64,1} = Array{Float64}(undef, numMP)

    num_surface_vertices::Int64 = length(lookup_tables[2].x)
    radii_diff::Array{Float64,1} = Array{Float64}(undef, num_surface_vertices)
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKBulkCustom(material_domain, bulk_data, pressure_diff, expected_pressure, lookup_tables[1], vtk_bulk_cells, pvd_bulk,
            t_index, t_record, t, joinpath(vtk_output_folder, "bulk"))
        writeVTKSurfaceCustom(surface_domain, surface_data, radii_diff, expected_radius, lookup_tables[2], vtk_surf_cells, pvd_surf,
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
    K::Float64 = arguments["bulk_modulus"]
    g::Float64 = arguments["gravity"]
    ρ::Float64 = arguments["density"]
    v::Float64 = arguments["dynamic_viscosity"]

    @threads for material_point::MaterialPoint in material_domain
        m::Float64 = material_point.volume_init * ρ

        material_point.mass = m
        material_point.density = ρ
        material_point.parameter_1 = K
        material_point.parameter_2 = v
        material_point.force_external = Vec2{Float64}(0.0, -m * g)
        material_point.velocity = ZERO_VEC2

        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices,
            material_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function setSurfaceParameters!(surface_domain::Vector{SurfacePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    γ::Float64 = arguments["surface_tension"]

    @threads for surface_point::SurfacePoint in surface_domain
        # Symmetry Lines x = 0 and y = 0
        if getSumOfVertexInDimension(surface_point, 1) < EPSILON || getSumOfVertexInDimension(surface_point, 2) < EPSILON
            surface_point.surface_tension = 0.0
            # Droplet-Air Surface
        else
            surface_point.surface_tension = γ
        end

        surface_point.connected_grid_array_length = getAllConnectedGrid!(surface_point.connected_grid_indices,
            surface_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function extraSetupDetails(arguments::Dict{String,Any}, expected_radius::Float64, expected_pressure::Float64)::String
    return """
    Notes
    - Units: grams (g), millimeters (mm), milliseconds (ms). Pressure is in MPa, force is in N.
    - The expected radius of the circle is $(expected_radius).
    - The expected pressure of the bulk is a constant $(expected_pressure)
    - Surface tension is loaded linearly over $(arguments["t_load"])
    - Symmetry boundaries on x=0 and y=0
    """
end
