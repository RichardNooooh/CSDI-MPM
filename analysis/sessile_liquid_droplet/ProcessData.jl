using Roots

function calculateExpectedRadiiAndPressure(arguments::Dict{String, Any})::Tuple{Float64, Float64}
    K::Float64 = arguments["bulk_modulus"]
    γ::Float64 = arguments["surface_tension"]
    s::Float64 = 2.0 # side length of initial full mesh

    f(r) = -K/2 * (π*r^2/s^2 - s^2/(π*r^2)) - γ/r
    expected_radius = find_zero(f, (0.5, 1.0)) # by default converges to machine tolerance
    expected_pressure = γ / expected_radius

    return expected_radius, expected_pressure
end


function writeVTKBulkCustom(material_domain::Vector{MaterialPoint},
                                bulk_data::BulkVTKData, plot_pressure_diff::Array{Float64, 1}, expected_pressure::Float64,
                                bulk_vtk_lookup::VTKLookup,
                                vtk_bulk_cells::Vector{MeshCell}, pvd_bulk::WriteVTK.CollectionFile,
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numMP::Int64 = length(material_domain)
    @threads for mp_index = 1:numMP
        material_point::MaterialPoint = material_domain[mp_index]

        # record all plotting information into these arrays
        bulk_data.pressure[mp_index] = material_point.pressure
        bulk_data.volume[mp_index] = material_point.volume
        bulk_data.volumeVerts[mp_index] = getQ4Area(material_point.vertices)
        bulk_data.density[mp_index] = material_point.density
        bulk_data.jacobian[mp_index] = det(material_point.deformation_gradient)
        bulk_data.velocity[mp_index, 1] = material_point.velocity[1]
        bulk_data.velocity[mp_index, 2] = material_point.velocity[2]

        bulk_data.stress[mp_index, 1, 1] = material_point.stress[1, 1]
        bulk_data.stress[mp_index, 1, 2] = material_point.stress[1, 2]
        bulk_data.stress[mp_index, 2, 1] = material_point.stress[2, 1]
        bulk_data.stress[mp_index, 2, 2] = material_point.stress[2, 2]

        bulk_data.deformgrad[mp_index, 1, 1] = material_point.deformation_gradient[1, 1]
        bulk_data.deformgrad[mp_index, 1, 2] = material_point.deformation_gradient[1, 2]
        bulk_data.deformgrad[mp_index, 2, 1] = material_point.deformation_gradient[2, 1]
        bulk_data.deformgrad[mp_index, 2, 2] = material_point.deformation_gradient[2, 2]

        plot_pressure_diff[mp_index] = material_point.pressure - expected_pressure
    end

    fillVertexCoordinates!(bulk_vtk_lookup, material_domain)

    # create vtk file
    vtk_grid(joinpath(vtk_output_folder, "bulk_timestep_$(Int64(t_index / t_record))"), 
                    bulk_vtk_lookup.x, bulk_vtk_lookup.y, bulk_vtk_lookup.z, vtk_bulk_cells) do vtk
        vtk["Pressure", VTKCellData()] =             bulk_data.pressure[:]
        vtk["Pressure Error", VTKCellData()] =       plot_pressure_diff[:]
        vtk["Volume", VTKCellData()] =               bulk_data.volume[:]
        vtk["Volume From Vertices", VTKCellData()] = bulk_data.volumeVerts[:]
        vtk["Density", VTKCellData()] =              bulk_data.density[:]
        vtk["Jacobian", VTKCellData()] =             bulk_data.jacobian[:]
        vtk["Velocity", VTKCellData()] =             @views (bulk_data.velocity[:, 1], bulk_data.velocity[:, 2], bulk_data.zero_cell)
        vtk["Stress", VTKCellData()] =               @views (bulk_data.stress[:, 1, 1], bulk_data.stress[:, 2, 1], bulk_data.zero_cell,
                                                            bulk_data.stress[:, 1, 2], bulk_data.stress[:, 2, 2], bulk_data.zero_cell,
                                                            bulk_data.zero_cell, bulk_data.zero_cell, bulk_data.zero_cell)
        vtk["Deformation Gradient", VTKCellData()] = @views (bulk_data.deformgrad[:, 1, 1], bulk_data.deformgrad[:, 2, 1], bulk_data.zero_cell,
                                                            bulk_data.deformgrad[:, 1, 2], bulk_data.deformgrad[:, 2, 2], bulk_data.zero_cell,
                                                            bulk_data.zero_cell, bulk_data.zero_cell, bulk_data.zero_cell)
        pvd_bulk[t] = vtk
    end

    return nothing
end


function writeVTKSurfaceCustom(surface_domain::Vector{SurfacePoint},
                                surface_data::SurfaceVTKData, plot_radii_diff::Array{Float64, 1}, expected_radius::Float64,
                                surface_vtk_lookup::VTKLookup,
                                vtk_surf_cells::Vector{MeshCell}, pvd_surf::WriteVTK.CollectionFile, 
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numSP::Int64 = length(surface_domain)
    @threads for sp_index = 1:numSP
        surface_point::SurfacePoint = surface_domain[sp_index]

        surface_data.surface_tension[sp_index] = surface_point.surface_tension
        surface_data.area[sp_index] = getL2Area(surface_point.vertices)
        
        normal_from_vertices = getL2Normal(surface_point.vertices)
        surface_data.normal[sp_index, 1] = normal_from_vertices[1]
        surface_data.normal[sp_index, 2] = normal_from_vertices[2]
    end

    fillVertexCoordinates!(surface_vtk_lookup, surface_domain)

    @threads for i in eachindex(surface_vtk_lookup.x)
        r = hypot(surface_vtk_lookup.x[i], surface_vtk_lookup.y[i])
        plot_radii_diff[i] = r - expected_radius
    end

    vtk_grid(joinpath(vtk_output_folder, "surf_timestep_$(Int64(t_index / t_record))"), 
                    surface_vtk_lookup.x, surface_vtk_lookup.y, surface_vtk_lookup.z, vtk_surf_cells) do vtk
        vtk["Surface Tension", VTKCellData()]      = surface_data.surface_tension[:]
        vtk["Area", VTKCellData()]                 = surface_data.area[:]
        vtk["Normal", VTKCellData()]               = @views (surface_data.normal[:, 1], surface_data.normal[:, 2], surface_data.zero_cell)
        vtk["Radius Error", VTKPointData()]        = plot_radii_diff[:]
        pvd_surf[t] = vtk
    end

    return nothing
end


function collect_surf_plotting_radiidiff!(surface_domain::Vector{SurfacePoint}, active_surface_domain::Vector{SurfacePoint},
                                    plot_surface_tension::Array{Float64, 1}, 
                                    plot_surf_cenradii::Array{Float64, 1}, plot_surf_cenradiidiff::Array{Float64, 1},
                                    plot_surface_area::Array{Float64, 1}, plot_surface_norm::Array{Float64, 2},
                                    vtk_surf_cells, pvd_surf, 
                                    t_index::Int64, t_record::Int64, t::Float64, expected_radius::Float64,
                                    vtk_output_folder::String)::Nothing
    vertex_dict::Dict{Int64, Vec2{Float64}} = Dict() # {VTK_Index, position}
    numSP::Int64 = length(surface_domain)
    @threads for sp_index = 1:numSP
        surface_point::SurfacePoint = surface_domain[sp_index]

        plot_surface_tension[sp_index] = surface_point.surface_tension

        if surface_point in active_surface_domain
            centroid::Vec2{Float64} = sum(surface_point.vertices) / 2
            centroid_radius::Float64 = norm(centroid)
            plot_surf_cenradii[sp_index] = centroid_radius
            plot_surf_cenradiidiff[sp_index] = centroid_radius - expected_radius
        else
            plot_surf_cenradii[sp_index] = 0.0
            plot_surf_cenradiidiff[sp_index] = 0.0
        end

        plot_surface_area[sp_index] = getL2Area(surface_point.vertices)
        normal_from_vertices = getL2Normal(surface_point.vertices)
        plot_surface_norm[sp_index, 1] = normal_from_vertices[1]
        plot_surface_norm[sp_index, 2] = normal_from_vertices[2]
    end

    for surface_point::SurfacePoint in surface_domain
        # find all unique vertex positions at this time step
        # surface_point.node_connectivity_vtk ensures the order is correct
        for (vertex_index, vtk_index) in enumerate(surface_point.node_connectivity_vtk)
            if haskey(vertex_dict, vtk_index)
                @assert vertex_dict[vtk_index] == surface_point.vertices[vertex_index]
            else
                # no lock needed since this is not multithreaded
                vertex_dict[vtk_index] = surface_point.vertices[vertex_index]
            end
        end
    end

    # create two lists for the x and y position of each vertex
    ordered_vtk_indices::Vector{Int64} = sort!(collect(keys(vertex_dict)))
    vtk_points_x::Vector{Float64} = Vector{Float64}(undef, length(ordered_vtk_indices))
    vtk_points_y::Vector{Float64} = Vector{Float64}(undef, length(ordered_vtk_indices))
    vtk_points_z::Vector{Float64} = zeros(length(ordered_vtk_indices))
    @threads for (list_index::Int64, vtk_key::Int64) in Folds.collect(enumerate(ordered_vtk_indices))
        vertex_position::Vec2{Float64} = vertex_dict[vtk_key]
        vtk_points_x[list_index] = vertex_position[1]
        vtk_points_y[list_index] = vertex_position[2]
    end

    # create vtk file
    zero_sp_array::Vector{Float64} = zeros(numSP)
    vtk_grid(string(vtk_output_folder, "surf_timestep_$(Int64(t_index / t_record))"), 
                    vtk_points_x, vtk_points_y, vtk_points_z, vtk_surf_cells) do vtk
        vtk["Surface Tension", VTKCellData()]      = plot_surface_tension[:]
        vtk["Radius", VTKCellData()]       = plot_surf_cenradii[:]
        vtk["RadiusDiff", VTKCellData()]   = plot_surf_cenradiidiff[:]
        vtk["Area", VTKCellData()]                 = plot_surface_area[:]
        vtk["Normal", VTKCellData()]               = @views (plot_surface_norm[:, 1], plot_surface_norm[:, 2], zero_sp_array)
        pvd_surf[t] = vtk
    end

    return nothing
end
