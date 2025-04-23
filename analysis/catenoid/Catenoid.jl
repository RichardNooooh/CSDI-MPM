include("../../src/3DMembrane/Driver.jl")
include("ProcessData.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments)
    @threads for grid_point::GridPoint in grid.points
        # top/bottom boundary
        if grid_point.position[3] < EPSILON || abs(grid_point.position[3] - 2.0) < EPSILON
            grid_point.is_fixed = Vec3{Bool}(true, true, true)
        end
    end

    membrane_domain::Vector{MembranePoint},
    membrane_lookup::VTKLookup = initializeBodyMesh(arguments, grid, setMembraneParameters!)

    # * ======================================
    # * Output Directory Setup
    # * ======================================
    destination_folder::String,
    vtk_output_folder::String,
    analysis_output_folder::String = initializeOutputDirectory(arguments)

    expected_catenoid_A::Float64 = calculatedExpectedCatenoidConstant()

    recordInitialSetupInfo(membrane_domain, grid,
        arguments, destination_folder, extraSetupDetails(arguments, expected_catenoid_A))

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    membrane_data::MembraneVTKData, pvd_membrane::WriteVTK.CollectionFile,
    vtk_membrane_cells::Vector{MeshCell} = initializeVTKMembraneData(membrane_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)
    writeVTIGrid(grid, destination_folder)

    t_record::Int64 = arguments["t_record"]
    num_membrane_vertices::Int64 = length(membrane_lookup.x)
    radii_diff::Array{Float64,1} = Array{Float64}(undef, num_membrane_vertices)
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKMembraneCustom(membrane_domain, membrane_data, radii_diff, expected_catenoid_A, membrane_lookup,
            vtk_membrane_cells, pvd_membrane,
            t_index, t_record, t, joinpath(vtk_output_folder, "membrane"))
        writeVTKGrid(grid, grid_data, pvd_grid,
            t_index, t_record, t, joinpath(vtk_output_folder, "grid"))
    end

    # * ======================================
    # * Setup "Active" Grid
    # * ======================================
    active_grid_points::Vector{GridPoint} = Folds.collect(grid_point
                                                          for grid_point::GridPoint in grid.points
                                                          if isActiveGridPoint(grid_point, arguments))
    @threads for gp_idx::Int64 in 1:length(grid.points) # mark these
        grid_data.active[gp_idx] = isActiveGridPoint(grid.points[gp_idx], arguments) ? 1.0 : 0.0
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
    compute_time::Int64 = runMPM!(arguments, grid, active_grid_points,
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

    @info("Recording All Positions Surface Point-Wise")
    open(analysis_output_folder * "vertex_list.tsv", "w") do file
        write(file, "x1\ty1\tz1\tx2\ty2\tz2\tx3\ty3\tz3\n")
        for membrane_point::MembranePoint in membrane_domain
            p1::Vec3{Float64}, p2::Vec3{Float64}, p3::Vec3{Float64} = membrane_point.vertices
            write(file, "$(p1[1])\t$(p1[2])\t$(p1[3])\t")
            write(file, "$(p2[1])\t$(p2[2])\t$(p2[3])\t")
            write(file, "$(p3[1])\t$(p3[2])\t$(p3[3])\n")
        end
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


function extraSetupDetails(arguments::Dict{String,Any}, expected_catenoid_A::Float64)::String
    return """
    Notes
    - Units: kilograms (kg), meters (m), seconds (s)
    - The expected constant (A) of the catenoid is $(expected_catenoid_A)
    - Surface tension is loaded linearly over $(arguments["t_load"])
    - Active grid nodes are points within r=[$(arguments["grid_r_min_active"]), $(arguments["grid_r_max_active"])] in cylindrical coordinates. 
    - Symmetry boundaries on x=0 and y=0 planes
    """
end

