include("../../src/3DMembrane/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments; num_extra_cells=2)

    membrane_domain::Vector{MembranePoint},
    membrane_lookup::VTKLookup = initializeBodyMesh(arguments, grid, setMembraneParameters!)

    rigid_domain::Vector{RigidPoint},
    rigid_lookup::VTKLookup = initializeRigidMesh(arguments, grid, setRigidBodyParameters!)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)
    if !isdir(joinpath(vtk_output_folder, "rigid"))
        mkdir(joinpath(vtk_output_folder, "rigid"))
    end

    recordInitialSetupInfo(membrane_domain, grid,
        arguments, destination_folder, extraSetupDetails(arguments))

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    membrane_data::MembraneVTKData, pvd_membrane::WriteVTK.CollectionFile,
    vtk_membrane_cells::Vector{MeshCell} = initializeVTKMembraneData(membrane_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)
    pvd_rigid::WriteVTK.CollectionFile,
    vtk_rigid_cells::Vector{MeshCell} = initializeVTKRigidData(rigid_domain, destination_folder)

    writeVTIGrid(grid, destination_folder)

    t_record::Int64 = arguments["t_record"]
    num_membrane_vertices::Int64 = length(membrane_lookup.x)
    radii::Array{Float64,1} = Array{Float64}(undef, num_membrane_vertices)
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKMembraneCustom(membrane_domain, membrane_data, radii, membrane_lookup,
            vtk_membrane_cells, pvd_membrane,
            t_index, t_record, t, joinpath(vtk_output_folder, "membrane"))
        writeVTKGrid(grid, grid_data, pvd_grid,
            t_index, t_record, t, joinpath(vtk_output_folder, "grid"))
        writeVTKRigid(rigid_domain, rigid_lookup, vtk_rigid_cells, pvd_rigid,
            t_index, t_record, t, joinpath(vtk_output_folder, "rigid"))
    end

    # * ======================================
    # * Setup "Active" Grid
    # * ======================================
    active_grid_points::Vector{GridPoint} = Folds.collect(grid_point
                                                          for grid_point::GridPoint in grid.points
                                                          if filterGridPoint(grid_point, arguments))
    @threads for gp_idx::Int64 in 1:length(grid.points) # mark these
        grid_data.active[gp_idx] = filterGridPoint(grid.points[gp_idx], arguments) ? 1.0 : 0.0
    end

    active_fixed_rigid_grid::Vector{GridPoint} = getActiveFixedGrid(grid.points, active_grid_points, rigid_domain)

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
    compute_time::Int64 = runMPM!(arguments, grid, active_grid_points, active_fixed_rigid_grid,
        membrane_domain,
        timeScaleFunction, recordVtkFiles!)

    # * ======================================
    # * Final Output Files
    # * ======================================
    @info("Saving PVD files")
    vtk_save(pvd_membrane)
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

    @info("Finished simulation! :)")
    open(analysis_output_folder * "COMPLETE", "w") do file
        write(file, "This is the receipt that the simulation is finished.\n")
    end
    return nothing
end


function setMembraneParameters!(membrane_domain::Vector{MembranePoint}, grid::Grid, arguments::Dict{String,Any})::Nothing
    γ::Float64 = arguments["surface_tension"]
    ρ::Float64 = arguments["area_density"]

    @threads for membrane_point::MembranePoint in membrane_domain
        m::Float64 = membrane_point.area * ρ
        membrane_point.mass = m
        membrane_point.surface_tension = γ

        membrane_point.connected_grid_array_length = getAllConnectedGrid!(membrane_point.connected_grid_indices,
            membrane_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end

function setRigidBodyParameters!(rigid_domain::Vector{RigidPoint}, grid::Grid)::Nothing
    @threads for rigid_point::RigidPoint in rigid_domain
        rigid_point.connected_grid_array_length = getAllConnectedGrid!(rigid_point.connected_grid_indices,
            rigid_point.vertices,
            grid.cell_length, grid.num_nodes)
    end
end

function extraSetupDetails(arguments::Dict{String, Any})::String
    return """
    Notes
    - Units: kilograms (kg), meters (m), seconds (s)
    - Rigid body mesh file: \"$(arguments["rigid_mesh_file"])\"
    - Surface tension is loaded linearly over $(arguments["t_load"])
    - Active grid nodes are points within a spherical radius of $(arguments["grid_sph_r_active"]) and cylindrical radii of $(arguments["grid_cyl_r_active"]).
    - Symmetry boundaries on x=0, y=0, and z=0 planes
    """
end
