function filterGridPoint(grid_point::GridPoint, arguments::Dict{String, Any})::Bool
    pos::Vec3{Float64} = grid_point.position
    cyl_r::Float64 = arguments["grid_cyl_r_active"]
    sph_r::Float64 = arguments["grid_sph_r_active"]
    if sqrt(pos[1]^2 + pos[2]^2) < cyl_r
        return false
    end
    if sqrt(pos[2]^2 + pos[3]^2) < cyl_r
        return false
    end
    if sqrt(pos[1]^2 + pos[3]^2) < cyl_r
        return false
    end

    if sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2) < sph_r
        return false
    end

    return true
end

function getActiveFixedGrid(grid_points::Vector{GridPoint}, active_grid_points::Vector{GridPoint},
    rigid_domain::Vector{RigidPoint})::Vector{GridPoint}

    numGP::Int64 = length(grid_points)
    marked_grid_points::Vector{Int64} = zeros(Int64, numGP)

    for rigid_point::RigidPoint in rigid_domain
        for connected_grid_index::Int8 in 1:rigid_point.connected_grid_array_length
            adj_grid_index::Int64 = rigid_point.connected_grid_indices[connected_grid_index]
            marked_grid_points[adj_grid_index] = 1
        end
    end

    rigid_grid_points::Vector{GridPoint} = Vector{GridPoint}(undef, 0)
    for i in eachindex(marked_grid_points)
        if marked_grid_points[i] == 1
            push!(rigid_grid_points, grid_points[i])
        end
    end

    active_rigid_grid_points = [rgp for rgp in rigid_grid_points if any(rgp === gp for gp in active_grid_points)]
    return active_rigid_grid_points
end


function writeVTKMembraneCustom(membrane_domain::Vector{MembranePoint},
                                membrane_data::MembraneVTKData, radii::Array{Float64, 1},
                                membrane_vtk_lookup::VTKLookup,
                                vtk_membrane_cells::Vector{MeshCell}, pvd_membrane::WriteVTK.CollectionFile, 
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numMP::Int64 = length(membrane_domain)
    @threads for mp_index = 1:numMP
        membrane_point::MembranePoint = membrane_domain[mp_index]

        membrane_data.surface_tension[mp_index] = membrane_point.surface_tension
        membrane_data.area[mp_index] = membrane_point.area
        
        normal::Vec3{Float64} = membrane_point.normal
        membrane_data.normal[mp_index, 1] = normal[1]
        membrane_data.normal[mp_index, 2] = normal[2]
        membrane_data.normal[mp_index, 3] = normal[3]

        velocity::Vec3{Float64} = membrane_point.velocity
        membrane_data.velocity[mp_index, 1] = velocity[1]
        membrane_data.velocity[mp_index, 2] = velocity[2]
        membrane_data.velocity[mp_index, 3] = velocity[3]
    end

    fillVertexCoordinates!(membrane_vtk_lookup, membrane_domain)

    @threads for i in eachindex(membrane_vtk_lookup.x)
        r = sqrt(membrane_vtk_lookup.x[i]^2 + membrane_vtk_lookup.y[i]^2 + membrane_vtk_lookup.z[i]^2)
        radii[i] = r
    end

    vtk_grid(joinpath(vtk_output_folder, "membrane_timestep_$(Int64(t_index / t_record))"), 
                    membrane_vtk_lookup.x, membrane_vtk_lookup.y, membrane_vtk_lookup.z, vtk_membrane_cells) do vtk
        vtk["Surface Tension", VTKCellData()] = membrane_data.surface_tension[:]
        vtk["Area", VTKCellData()]            = membrane_data.area[:]
        vtk["Normal", VTKCellData()]          = @views (membrane_data.normal[:, 1], membrane_data.normal[:, 2], membrane_data.normal[:, 3])
        vtk["Velocity", VTKCellData()]        = @views (membrane_data.velocity[:, 1], membrane_data.velocity[:, 2], membrane_data.velocity[:, 3])
        vtk["Radius", VTKPointData()]         = radii[:]
        pvd_membrane[t] = vtk
    end

    return nothing
end