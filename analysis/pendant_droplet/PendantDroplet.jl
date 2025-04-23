include("../../src/3D/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments; num_extra_cells=2)

    for grid_point::GridPoint in grid.points
        # ceiling for droplet
        if abs(grid_point.position[3] - 3.0) < EPSILON
            grid_point.is_fixed = Vec3{Bool}(true, true, true)
        end
    end

    material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    active_surface_domain::Vector{SurfacePoint},
    lookup_tables::Tuple{VTKLookup,VTKLookup} = initializeBodyMeshes(arguments, grid, setBulkParameters!, setSurfaceParameters!; get_active_subset=false)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)

    recordInitialSetupInfo(material_domain, surface_domain, active_surface_domain,
        grid, arguments, destination_folder, extraSetupDetails(arguments))

    writeVTIGrid(grid, destination_folder)

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    bulk_data::BulkVTKData, pvd_bulk::WriteVTK.CollectionFile,
    vtk_bulk_cells::Vector{MeshCell} = initializeVTKBulkData(material_domain, destination_folder)
    surface_data::SurfaceVTKData, pvd_surf::WriteVTK.CollectionFile,
    vtk_surf_cells::Vector{MeshCell} = initializeVTKSurfaceData(surface_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)

    t_record::Int64 = arguments["t_record"]
    num_surface_vertices::Int64 = length(lookup_tables[2].x)
    radii::Array{Float64,1} = Array{Float64}(undef, num_surface_vertices)
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKBulk(material_domain, bulk_data, lookup_tables[1], vtk_bulk_cells, pvd_bulk,
            t_index, t_record, t, joinpath(vtk_output_folder, "bulk"))
        writeVTKSurfaceCustom(surface_domain, surface_data, radii, lookup_tables[2], vtk_surf_cells, pvd_surf,
            t_index, t_record, t, joinpath(vtk_output_folder, "surface"))
        writeVTKGrid(grid, grid_data, pvd_grid,
            t_index, t_record, t, joinpath(vtk_output_folder, "grid"))
    end

    # * ======================================
    # * Time Scale Functions
    # * ======================================
    t_load_g_start::Float64, t_load_g_max::Float64 = arguments["t_load_g_start"], arguments["t_load_g_max"]
    t_load_st_max::Float64, t_load_st_decay_start::Float64, t_load_st_decay_factor::Float64 =
        arguments["t_load_st_max"], arguments["t_load_st_decay_start"], arguments["t_load_st_decay_factor"]

    function timeScaleGravity(t::Float64)::Float64
        return maximum((0.0, minimum((t_load_g_max, t - t_load_g_start)) / t_load_g_max))
    end

    function timeScaleSurfaceTension(t::Float64)::Float64
        if t < t_load_st_max
            return t / t_load_st_max
        elseif t < t_load_st_decay_start
            return 1.0
        else
            k::Float64 = t_load_st_decay_factor
            return exp(-k * (t - t_load_st_decay_start))
        end
    end

    # * ======================================
    # * Running MPM Algorithm (Steady State Stage)
    # * ======================================
    dt::Float64, t_final::Float64, v_alpha::Float64 = arguments["t_delta"], arguments["t_final"], arguments["v_alpha"]

    arguments["t_final"] = arguments["t_load_st_decay_start"]
    arguments["v_alpha"] = arguments["v_alpha_initial"]
    @info("Stage 1: Compressive loading to steady state")
    stage_1_time::Int64 = runMPM!(arguments, grid,
        material_domain, surface_domain, active_surface_domain,
        timeScaleGravity, timeScaleSurfaceTension, recordVtkFiles!)

    # * ======================================
    # * Running MPM Algorithm (Decay Stage)
    # * ======================================
    stage_2_start_index = Int64(arguments["t_load_st_decay_start"] / dt) + 1
    arguments["t_init"] = arguments["t_load_st_decay_start"] + dt
    arguments["t_final"] = t_final
    arguments["v_alpha"] = v_alpha
    setDynVisc(material_domain, arguments)

    @info("Stage 2: Surface Tension Exponential Decay")
    stage_2_time::Int64 = runMPM!(arguments, grid,
        material_domain, surface_domain, active_surface_domain,
        timeScaleGravity, timeScaleSurfaceTension, recordVtkFiles!;
        start_index=stage_2_start_index)

    compute_time::Int64 = stage_1_time + stage_2_time

    # * ======================================
    # * Final Output Files
    # * ======================================
    @info("Saving PVD files")
    vtk_save(pvd_bulk)
    vtk_save(pvd_surf)
    vtk_save(pvd_grid)

    @info("Recording compute time")
    open(joinpath(analysis_output_folder, "compute_time_ms.txt"), "w") do file
        write(file, "$compute_time\n")
    end

    @info("Finished simulation! :)")
    open(joinpath(analysis_output_folder, "COMPLETE"), "w") do file
        write(file, "This is the receipt that the simulation is finished.\n")
    end

    return nothing
end


function setDynVisc(material_domain::Vector{MaterialPoint}, arguments::Dict{String, Any})::Nothing
    µ::Float64 = arguments["dynamic_viscosity"]
    @threads for material_point::MaterialPoint in material_domain
        material_point.parameter_2 = µ
    end
end


function setBulkParameters!(material_domain::Vector{MaterialPoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    K::Float64 = arguments["bulk_modulus"]
    g::Float64 = arguments["gravity"]
    ρ::Float64 = arguments["density"]
    translate::Vector{Float64} = arguments["translation_vector"]

    @threads for material_point::MaterialPoint in material_domain
        m::Float64 = material_point.volume_init * ρ

        for vertex_index::Int64 in eachindex(material_point.vertices)
            material_point.vertices[vertex_index] += Vec3{Float64}(translate[1], translate[2], translate[3])
        end

        material_point.mass = m
        material_point.density = ρ
        material_point.parameter_1 = K
        material_point.parameter_2 = 0.0
        material_point.force_external = Vec3{Float64}(0.0, 0.0, -m * g)
        material_point.velocity = ZERO_VEC3

        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices,
            material_point.vertices,
            grid.cell_length, grid.num_nodes)
    end

end


function setSurfaceParameters!(surface_domain::Vector{SurfacePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    γ::Float64 = arguments["surface_tension"]
    translate::Vector{Float64} = arguments["translation_vector"]

    @threads for surface_point::SurfacePoint in surface_domain
        # Symmetry Planes x = 0, y = 0 and Ceiling z = 1
        if getSumOfVertexInDimension(surface_point, 1) < EPSILON || getSumOfVertexInDimension(surface_point, 2) < EPSILON || abs(getSumOfVertexInDimension(surface_point, 3) / 3.0 - 3.0) < EPSILON
            surface_point.surface_tension = 0.0
        else
            surface_point.surface_tension = γ
        end

        for vertex_index::Int64 in eachindex(surface_point.vertices)
            surface_point.vertices[vertex_index] += Vec3{Float64}(translate[1], translate[2], translate[3])
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
    - Translates mesh by $(arguments["translation_vector"]).
    - Time Scaling Parameters for Gravity:
      - Time to start:         $(arguments["t_load_g_start"])
      - Time to reach maximum: $(arguments["t_load_g_start"] + arguments["t_load_g_max"])
    - Time Scaling Parameters for Surface Tension:
      - Time to reach maximum: $(arguments["t_load_st_max"])
      - Time to start decay:   $(arguments["t_load_st_decay_start"])
      - Exponential factor:    $(arguments["t_load_st_decay_factor"])
    - v_alpha above is for the 2nd exp. decay stage. Initial v_alpha is $(arguments["v_alpha_initial"])
    """
end
