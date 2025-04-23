
function writeVTKSurfaceCustom(surface_domain::Vector{SurfacePoint},
                                surface_data::SurfaceVTKData, radii::Array{Float64, 1}, surface_vtk_lookup::VTKLookup,
                                vtk_surf_cells::Vector{MeshCell}, pvd_surf::WriteVTK.CollectionFile, 
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numSP::Int64 = length(surface_domain)
    @threads for sp_index = 1:numSP
        surface_point::SurfacePoint = surface_domain[sp_index]

        surface_data.surface_tension[sp_index] = surface_point.surface_tension
        surface_data.area[sp_index] = getTri3Area(surface_point.vertices)
        
        normal_from_vertices = getTri3Normal(surface_point.vertices)
        surface_data.normal[sp_index, 1] = normal_from_vertices[1]
        surface_data.normal[sp_index, 2] = normal_from_vertices[2]
        surface_data.normal[sp_index, 3] = normal_from_vertices[3]
    end

    fillVertexCoordinates!(surface_vtk_lookup, surface_domain)

    @threads for i in eachindex(surface_vtk_lookup.x)
        r = sqrt(surface_vtk_lookup.x[i]^2 + surface_vtk_lookup.y[i]^2)
        radii[i] = r
    end

    vtk_grid(joinpath(vtk_output_folder, "surf_timestep_$(Int64(t_index / t_record))"), 
                    surface_vtk_lookup.x, surface_vtk_lookup.y, surface_vtk_lookup.z, vtk_surf_cells) do vtk
        vtk["Surface Tension", VTKCellData()]      = surface_data.surface_tension[:]
        vtk["Area", VTKCellData()]                 = surface_data.area[:]
        vtk["Normal", VTKCellData()]               = @views (surface_data.normal[:, 1], surface_data.normal[:, 2], surface_data.normal[:, 3])
        vtk["Radius", VTKPointData()]              = radii[:]
        pvd_surf[t] = vtk
    end

    return nothing
end
