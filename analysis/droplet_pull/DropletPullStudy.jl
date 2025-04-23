include("../../src/2D/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments)

    material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint}, UNUSED::Vector{SurfacePoint},
    lookup_tables::Tuple{VTKLookup,VTKLookup} = initializeBodyMeshes(arguments, grid, 
        setBulkParameters!, setSurfaceParameters!; get_active_subset=false)

    # rigid body for the 2nd stage
    rigid_domain::Vector{RigidPoint}, 
    rigid_lookup::VTKLookup = initializeRigidMesh(arguments, grid, setRigidParameters!)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)
    if !isdir(joinpath(vtk_output_folder, "rigid"))
        mkdir(joinpath(vtk_output_folder, "rigid"))
    end

    recordInitialSetupInfo(material_domain, surface_domain, surface_domain,
        grid, arguments, destination_folder, extraSetupDetails(arguments))

    # * ======================================
    # * Time Scale Functions
    # * ======================================
    t_load::Float64 = arguments["t_load"]
    function timeScaleFunction(t::Float64)::Float64
        return minimum((t_load, t)) / t_load
    end

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    bulk_data::BulkVTKData, pvd_bulk::WriteVTK.CollectionFile,
    vtk_bulk_cells::Vector{MeshCell} = initializeVTKBulkData(material_domain, destination_folder)
    surface_data::SurfaceVTKData, pvd_surf::WriteVTK.CollectionFile,
    vtk_surf_cells::Vector{MeshCell} = initializeVTKSurfaceData(surface_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)
    pvd_rigid::WriteVTK.CollectionFile,
    vtk_rigid_cells::Vector{MeshCell} = initializeVTKRigidData(rigid_domain, destination_folder)
    writeVTIGrid(grid, destination_folder)

    # get grid points on axis of symmetry
    symAxisGridPoints::Dict{UInt32,Tuple{GridPoint,GridPoint}} = getTestGridpointsOnSymmetry(grid, 2.75)
    
    left_material_domain::Vector{MaterialPoint}, 
    right_material_domain::Vector{MaterialPoint},
    left_surface_domain::Vector{SurfacePoint},
    right_surface_domain::Vector{SurfacePoint} = getLeftRightDroplets(material_domain, surface_domain, 2.75)

    numGP::Int64 = length(grid.points)
    plot_left_force::Array{Float64, 2} = Array{Float64}(undef, numGP, 2);
    plot_right_force::Array{Float64, 2} = Array{Float64}(undef, numGP, 2);

    plot_force_over_time::Vector{Tuple{Float64, Vec2{Float64}, Vec2{Float64}}} = Vector{Tuple{Float64, Vec2{Float64}, Vec2{Float64}}}(undef, 0);

    t_record::Int64 = arguments["t_record"]
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKBulk(material_domain, bulk_data, lookup_tables[1], vtk_bulk_cells, pvd_bulk,
            t_index, t_record, t, joinpath(vtk_output_folder, "bulk"))
        writeVTKSurface(surface_domain, surface_data, lookup_tables[2], vtk_surf_cells, pvd_surf,
            t_index, t_record, t, joinpath(vtk_output_folder, "surface"))

        # reset grid
        for (gp_left::GridPoint, gp_right::GridPoint) in values(symAxisGridPoints)
            gp_left.mass = 0.0
            gp_left.momentum = ZERO_VEC2
            gp_left.force = ZERO_VEC2

            gp_right.mass = 0.0
            gp_right.momentum = ZERO_VEC2
            gp_right.force = ZERO_VEC2
        end
        # do material/surface to grid, one side at a time then record
        t_scale::Float64 = timeScaleFunction(t)
        oneSidedGridForceCalculation!(grid, false, left_material_domain, left_surface_domain, symAxisGridPoints, t_scale)
        oneSidedGridForceCalculation!(grid, true, right_material_domain, right_surface_domain, symAxisGridPoints, t_scale)
        writeVTKGridCustom!(grid, grid_data, symAxisGridPoints,
            plot_left_force, plot_right_force, plot_force_over_time,
            pvd_grid, t_index, t_record, t, joinpath(vtk_output_folder, "grid"))

        writeVTKRigidBulk(rigid_domain, rigid_lookup, vtk_rigid_cells, pvd_rigid,
            t_index, t_record, t, joinpath(vtk_output_folder, "rigid"))
    end

    # * ======================================
    # * Running MPM Algorithm (Steady State Stage)
    # * ======================================
    dt::Float64, t_final::Float64, v_alpha::Float64 = arguments["t_delta"], arguments["t_final"], arguments["v_alpha"]

    arguments["t_final"] = arguments["t_velocity_start"]
    arguments["v_alpha"] = arguments["v_alpha_initial"]
    @info("Stage 1: Compressive loading to steady state")
    stage_1_time::Int64 = runMPM!(arguments, grid,
        material_domain, surface_domain, surface_domain,
        timeScaleFunction, timeScaleFunction, recordVtkFiles!)

    # * ======================================
    # * Running MPM Algorithm (Pull Stage)
    # * ======================================
    stage_2_start_index = Int64(arguments["t_velocity_start"] / dt) + 1
    arguments["t_init"] = arguments["t_velocity_start"] + dt
    arguments["t_final"] = t_final
    arguments["v_alpha"] = v_alpha
    @info("Stage 2: Pull apart droplets")
    stage_2_time::Int64 = runMPM!(arguments, grid,
        material_domain, surface_domain, surface_domain,
        rigid_domain,
        timeScaleFunction, timeScaleFunction, recordVtkFiles!;
        start_index=stage_2_start_index)

    compute_time::Int64 = stage_1_time + stage_2_time

    # * ======================================
    # * Final Output Files
    # * ======================================
    @info("Saving PVD files")
    vtk_save(pvd_bulk)
    vtk_save(pvd_surf)
    vtk_save(pvd_grid)
    vtk_save(pvd_rigid)

    @info("Recording compute time")
    open(joinpath(analysis_output_folder, "compute_time_ms.txt"), "w") do file
        write(file, "$compute_time\n")
    end

    @info("Recording the left/right total forces on left and right droplets along symmetry axis.")
    open(joinpath(analysis_output_folder, "symmetry_forces_time.tsv"), "w") do file
        write(file, "time\tleft_force(x)\tleft_force(y)\tright_force(x)\tright_force(y)\n")
        for (t::Float64, left_force::Vec2{Float64}, right_force::Vec2{Float64}) in plot_force_over_time
            write(file, "$(t)\t$(left_force[1])\t$(left_force[2])\t$(right_force[1])\t$(right_force[2])\n")
        end
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


function setRigidParameters!(rigid_domain::Vector{RigidPoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    v::Float64 = arguments["velocity"]

    translate1::Vector{Float64} = arguments["translation_vector_1"]
    translate2::Vector{Float64} = arguments["translation_vector_2"]
    @threads for rigid_point::RigidPoint in rigid_domain
        centroid::Vec2{Float64} = sum(rigid_point.vertices) / 4
        average_x::Float64 = centroid[1]
        if average_x < 3.25
            rigid_point.velocity = Vec2{Float64}(-v, 0.0)
            for vertex_index::Int64 in eachindex(rigid_point.vertices)
                rigid_point.vertices[vertex_index] += Vec2{Float64}(translate1[1], translate1[2])
            end
        else
            rigid_point.velocity = Vec2{Float64}(v, 0.0)
            for vertex_index::Int64 in eachindex(rigid_point.vertices)
                rigid_point.vertices[vertex_index] += Vec2{Float64}(translate2[1], translate2[2])
            end
        end

        rigid_point.connected_grid_array_length = getAllConnectedGrid!(rigid_point.connected_grid_indices,
            rigid_point.vertices,
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
    - Pull start time of droplets: $(arguments["t_velocity_start"])
    - Velocity vector magnitude: = $(v)
    - Rigid body mesh file: $(arguments["rigid_mesh_file"])
    - Droplet 1 translation vector from droplet center, (1.25, 1.25): $(arguments["translation_vector_1"])
    - Droplet 2 translation vector from droplet center, (4.25, 1.25): $(arguments["translation_vector_2"])
    - Droplet edges are $(1.0 - arguments["translation_vector_1"][1] + arguments["translation_vector_2"][1]) away from each other.
    """
end

