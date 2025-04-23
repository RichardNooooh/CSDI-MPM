include("../../src/2D/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments)

    material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    UNUSED::Vector{SurfacePoint},
    lookup_tables::Tuple{VTKLookup,VTKLookup} = initializeBodyMeshes(arguments, grid,
        setBulkParameters!, setSurfaceParameters!; get_active_subset=false)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)

    recordInitialSetupInfo(material_domain, surface_domain, surface_domain,
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
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKBulk(material_domain, bulk_data, lookup_tables[1], vtk_bulk_cells, pvd_bulk,
            t_index, t_record, t, joinpath(vtk_output_folder, "bulk"))
        writeVTKSurface(surface_domain, surface_data, lookup_tables[2], vtk_surf_cells, pvd_surf,
            t_index, t_record, t, joinpath(vtk_output_folder, "surface"))
        writeVTKGrid(grid, grid_data, pvd_grid,
            t_index, t_record, t, joinpath(vtk_output_folder, "grid"))
    end

    # * ======================================
    # * Time Scale Functions
    # * ======================================
    t_load::Float64 = arguments["t_load"]
    function timeScaleFunction(t::Float64)::Float64
        return minimum((t_load, t)) / t_load
    end

    # * ======================================
    # * Running MPM Algorithm (Steady State Stage)
    # * ======================================
    dt::Float64, t_final::Float64, v_alpha::Float64 = arguments["t_delta"], arguments["t_final"], arguments["v_alpha"]

    arguments["t_final"] = arguments["t_velocity_start"]
    arguments["v_alpha"] = arguments["v_alpha_initial"]
    @info("Stage 1: Compressive loading to steady state")
    stage_1_time::Int64 = runMPM!(arguments, grid, material_domain, surface_domain, surface_domain,
        timeScaleFunction, timeScaleFunction, recordVtkFiles!)

    # * ======================================
    # * Running MPM Algorithm (Launch Stage)
    # * ======================================
    setVelocity(material_domain, arguments)
    stage_2_start_index = Int64(arguments["t_velocity_start"] / dt) + 1
    arguments["t_init"] = arguments["t_velocity_start"] + dt
    arguments["t_final"] = t_final
    arguments["v_alpha"] = v_alpha
    @info("Stage 2: Launch to collision")
    stage_2_time::Int64 = runMPM!(arguments, grid, material_domain, surface_domain, surface_domain,
        timeScaleFunction, timeScaleFunction, recordVtkFiles!; start_index=stage_2_start_index)

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

function setVelocity(material_domain::Vector{MaterialPoint}, arguments::Dict{String,Any})::Nothing
    vel::Float64 = arguments["velocity"]

    @threads for material_point::MaterialPoint in material_domain
        centroid::Vec2{Float64} = sum(material_point.vertices) / 4
        average_x::Float64 = centroid[1]
        if average_x < 3.25
            material_point.velocity = Vec2{Float64}(vel, vel)
        else
            material_point.velocity = Vec2{Float64}(-vel, -vel)
        end
    end
end


function setBulkParameters!(material_domain::Vector{MaterialPoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    E::Float64 = arguments["elastic_modulus"]
    v::Float64 = arguments["poisson_ratio"]
    g::Float64 = arguments["gravity"]
    ρ::Float64 = arguments["density"]

    λ::Float64, μ::Float64 = getLaméParameters(E, v)

    translate1::Vector{Float64} = arguments["translation_vector_1"]
    translate2::Vector{Float64} = arguments["translation_vector_2"]
    @threads for material_point::MaterialPoint in material_domain
        m::Float64 = material_point.volume_init * ρ

        centroid::Vec2{Float64} = sum(material_point.vertices) / 4
        average_x::Float64 = centroid[1]
        if average_x < 3.25
            for vertex_index::Int64 in eachindex(material_point.vertices)
                material_point.vertices[vertex_index] += Vec2{Float64}(translate1[1], translate1[2])
            end
        else
            for vertex_index::Int64 in eachindex(material_point.vertices)
                material_point.vertices[vertex_index] += Vec2{Float64}(translate2[1], translate2[2])
            end
        end

        material_point.mass = m
        material_point.density = ρ
        material_point.force_external = Vec2{Float64}(0.0, -m * g)

        material_point.parameter_1 = λ
        material_point.parameter_2 = μ

        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices,
            material_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end


function setSurfaceParameters!(surface_domain::Vector{SurfacePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    γ::Float64 = arguments["surface_tension"]

    translate1::Vector{Float64} = arguments["translation_vector_1"]
    translate2::Vector{Float64} = arguments["translation_vector_2"]
    @threads for surface_point::SurfacePoint in surface_domain
        surface_point.surface_tension = γ

        centroid::Vec2{Float64} = sum(surface_point.vertices) / 2
        average_x::Float64 = centroid[1]

        if average_x < 3.25
            for vertex_index::Int64 in eachindex(surface_point.vertices)
                surface_point.vertices[vertex_index] += Vec2{Float64}(translate1[1], translate1[2])
            end
        else
            for vertex_index::Int64 in eachindex(surface_point.vertices)
                surface_point.vertices[vertex_index] += Vec2{Float64}(translate2[1], translate2[2])
            end
        end

        surface_point.connected_grid_array_length = getAllConnectedGrid!(surface_point.connected_grid_indices,
            surface_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end

function extraSetupDetails(arguments::Dict{String,Any})::String
    v::Float64 = arguments["velocity"]
    return """
    Notes
    - Units: grams (g), millimeters (mm), milliseconds (ms). Pressure is in MPa, force is in N.
    - Surface tension is loaded linearly over $(arguments["t_load"])
    - Velocity alpha (above) is for stage 2. Stage 1 has velocity alpha at $(arguments["v_alpha_initial"]).
    - Launch time of droplets: $(arguments["t_velocity_start"])
    - Velocity vector magnitude: sqrt($v^2 + $v^2) = $(sqrt(v^2+v^2))
    - Droplet 1 translation vector from droplet center, (1.25, 1.25): $(arguments["translation_vector_1"])
    - Droplet 2 translation vector from droplet center, (4.25, 1.25): $(arguments["translation_vector_2"])
    """
end

