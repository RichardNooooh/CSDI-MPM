include("../src/ValidateArgs.jl")

_arguments = Dict{String,Any}([
    ("simulation_name", "Fork Membrane"),
    ("mesh_file", "meshes/3D/forktubes_msh36.msh"),
    ("grid_dims", [4.5, 2.5, 1.0]),
    ("grid_num_cells", [180, 100, 40]),
    ("dest", "output/forktubes/"),
    ("output_name", "h40/"),
    ("t_init", 0.0),
    ("t_delta", 1e-3),
    ("t_final", 500.0),
    ("t_load", 100),
    ("t_record", 2000),
    ("v_alpha", 0.0),
    ("area_density", 1.0),
    ("surface_tension", 1.0)
])
validateMembraneArguments(_arguments)

include("../src/3DMembrane/Driver.jl")

function run(arguments::Dict{String,Any})::Nothing
    # * ======================================
    # * Grid and Point Domain Initialization
    # * ======================================
    grid::Grid = initializeGrid(arguments)
    @threads for grid_point::GridPoint in grid.points
        # top/bottom boundary
        if grid_point.position[3] < EPSILON || abs(grid_point.position[3] - 1.0) < EPSILON
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

    recordInitialSetupInfo(membrane_domain, grid,
        arguments, destination_folder, extraSetupDetails(arguments))

    # * ======================================
    # * VTK Output Setup
    # * ======================================
    membrane_data::MembraneVTKData, pvd_membrane::WriteVTK.CollectionFile,
    vtk_membrane_cells::Vector{MeshCell} = initializeVTKMembraneData(membrane_domain, destination_folder)
    grid_data::GridVTKData, pvd_grid::WriteVTK.CollectionFile = initializeVTKGridData(grid, destination_folder)
    writeVTIGrid(grid, destination_folder)

    t_record::Int64 = arguments["t_record"]
    function recordVtkFiles!(t::Float64, t_index::Int64)
        writeVTKMembrane(membrane_domain, membrane_data, membrane_lookup,
            vtk_membrane_cells, pvd_membrane,
            t_index, t_record, t, joinpath(vtk_output_folder, "membrane"))
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
    compute_time::Int64 = runMPM!(arguments, grid,
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



function extraSetupDetails(arguments::Dict{String,Any})::String
    return """
    Notes
    - Units: kilograms (kg), meters (m), seconds (s)
    - Surface tension is loaded linearly over $(arguments["t_load"])
    """
end



run(_arguments)
