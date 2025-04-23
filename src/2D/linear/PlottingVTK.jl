
function collect_bulk_plotting!(material_domain::Vector{MaterialPoint},
                                plot_pressure::Array{Float64, 1}, plot_volume::Array{Float64, 1},
                                plot_density::Array{Float64, 1},
                                plot_jacobian::Array{Float64, 1}, plot_vel::Array{Float64, 2},
                                plot_stress::Array{Float64, 3}, plot_strain::Array{Float64, 3},
                                plot_deformgrad::Array{Float64, 3},
                                vtk_bulk_cells, pvd_bulk,
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numMP::Int64 = length(material_domain)
    @threads for mp_index = 1:numMP
        material_point::MaterialPoint = material_domain[mp_index]

        # record all plotting information into these arrays
        plot_pressure[mp_index] = material_point.pressure
        plot_volume[mp_index] = material_point.volume
        # plot_volumeVerts[mp_index] = getQ4Area(material_point.vertices)
        plot_density[mp_index] = material_point.density
        plot_jacobian[mp_index] = det(material_point.deformation_gradient)
        plot_vel[mp_index, 1] = material_point.velocity[1]
        plot_vel[mp_index, 2] = material_point.velocity[2]

        plot_stress[mp_index, 1, 1] = material_point.stress[1, 1]
        plot_stress[mp_index, 1, 2] = material_point.stress[1, 2]
        plot_stress[mp_index, 2, 1] = material_point.stress[2, 1]
        plot_stress[mp_index, 2, 2] = material_point.stress[2, 2]

        plot_strain[mp_index, 1, 1] = material_point.strain[1, 1]
        plot_strain[mp_index, 1, 2] = material_point.strain[1, 2]
        plot_strain[mp_index, 2, 1] = material_point.strain[2, 1]
        plot_strain[mp_index, 2, 2] = material_point.strain[2, 2]

        plot_deformgrad[mp_index, 1, 1] = material_point.deformation_gradient[1, 1]
        plot_deformgrad[mp_index, 1, 2] = material_point.deformation_gradient[1, 2]
        plot_deformgrad[mp_index, 2, 1] = material_point.deformation_gradient[2, 1]
        plot_deformgrad[mp_index, 2, 2] = material_point.deformation_gradient[2, 2]
    end

    vtk_points_x::Vector{Float64} = Vector{Float64}(undef, numMP)
    vtk_points_y::Vector{Float64} = Vector{Float64}(undef, numMP)
    vtk_points_z::Vector{Float64} = zeros(numMP)
    @threads for mp_index::Int64 in eachindex(material_domain)
        position::Vec2{Float64} = material_domain[mp_index].position
        vtk_points_x[mp_index] = position[1]
        vtk_points_y[mp_index] = position[2]
    end

    # create vtk file
    zero_mp_array::Vector{Float64} = zeros(numMP)
    vtk_grid(string(vtk_output_folder, "bulk_timestep_$(Int64(t_index / t_record))"), 
                    vtk_points_x, vtk_points_y, vtk_points_z, vtk_bulk_cells) do vtk
        vtk["Pressure", VTKPointData()] =             plot_pressure[:]
        vtk["Volume", VTKPointData()] =               plot_volume[:]
        # vtk["Volume From Vertices", VTKPointData()] = plot_volumeVerts[:]
        vtk["Density", VTKPointData()] =              plot_density[:]
        vtk["Jacobian", VTKPointData()] =             plot_jacobian[:]
        vtk["Velocity", VTKPointData()] =             @views (plot_vel[:, 1], plot_vel[:, 2], zero_mp_array)
        vtk["Stress", VTKPointData()] =               @views (plot_stress[:, 1, 1], plot_stress[:, 2, 1], zero_mp_array,
                                                            plot_stress[:, 1, 2], plot_stress[:, 2, 2], zero_mp_array,
                                                            zero_mp_array, zero_mp_array, zero_mp_array)
        vtk["Strain", VTKPointData()] =               @views (plot_strain[:, 1, 1], plot_strain[:, 2, 1], zero_mp_array,
                                                            plot_strain[:, 1, 2], plot_strain[:, 2, 2], zero_mp_array,
                                                            zero_mp_array, zero_mp_array, zero_mp_array)
        vtk["Deformation Gradient", VTKPointData()] = @views (plot_deformgrad[:, 1, 1], plot_deformgrad[:, 2, 1], zero_mp_array,
                                                            plot_deformgrad[:, 1, 2], plot_deformgrad[:, 2, 2], zero_mp_array,
                                                            zero_mp_array, zero_mp_array, zero_mp_array)
        pvd_bulk[t] = vtk
    end

    return nothing
end

function collect_surf_plotting!(surface_domain::Vector{SurfacePoint},
                                plot_surface_tension::Array{Float64, 1}, plot_surface_area::Array{Float64, 1}, 
                                plot_surface_norm::Array{Float64, 2}, 
                                vtk_surf_cells, pvd_surf, 
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numSP::Int64 = length(surface_domain)
    @threads for sp_index = 1:numSP
        surface_point::SurfacePoint = surface_domain[sp_index]

        plot_surface_tension[sp_index] = surface_point.surface_tension
        plot_surface_area[sp_index] = surface_point.area
        plot_surface_norm[sp_index, 1] = surface_point.normal[1]
        plot_surface_norm[sp_index, 2] = surface_point.normal[2]
    end

    # create two lists for the x and y position of each vertex
    vtk_points_x::Vector{Float64} = Vector{Float64}(undef, numSP)
    vtk_points_y::Vector{Float64} = Vector{Float64}(undef, numSP)
    vtk_points_z::Vector{Float64} = zeros(numSP)
    @threads for sp_index::Int64 in eachindex(surface_domain)
        position::Vec2{Float64} = surface_domain[sp_index].position
        vtk_points_x[sp_index] = position[1]
        vtk_points_y[sp_index] = position[2]
    end

    # create vtk file
    zero_sp_array::Vector{Float64} = zeros(numSP)
    vtk_grid(string(vtk_output_folder, "surf_timestep_$(Int64(t_index / t_record))"), 
                    vtk_points_x, vtk_points_y, vtk_points_z, vtk_surf_cells) do vtk
        vtk["Surface Tension", VTKCellData()]      = plot_surface_tension[:]
        vtk["Area", VTKCellData()]                 = plot_surface_area[:]
        vtk["Normal", VTKCellData()]               = @views (plot_surface_norm[:, 1], plot_surface_norm[:, 2], zero_sp_array)
        pvd_surf[t] = vtk
    end

    return nothing
end

function collect_grid_plotting!(grid::Grid,
                                plot_gridmass::Array{Float64, 1},
                                plot_gridmomentum::Array{Float64, 2},
                                plot_gridforce::Array{Float64, 2},
                                pvd_grid,
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numGP::Int64 = length(grid.points)
    for gp_index = 1:numGP
        grid_point::GridPoint = grid.points[gp_index]

        # record all plotting information into these arrays
        plot_gridmass[gp_index] = grid_point.mass
        plot_gridmomentum[gp_index, 1] = grid_point.momentum[1]
        plot_gridmomentum[gp_index, 2] = grid_point.momentum[2]
        plot_gridforce[gp_index, 1] = grid_point.force[1]
        plot_gridforce[gp_index, 2] = grid_point.force[2]
    end

    # create vtk file
    zero_gp_array::Vector{Float64} = zeros(numGP)
    vtk_grid(string(vtk_output_folder, "grid_timestep_$(Int64(t_index / t_record))"), 
             0.0:grid.cell_length[1]:grid.grid_length[1],
             0.0:grid.cell_length[2]:grid.grid_length[2]) do vtk
        vtk["Mass"] = plot_gridmass[:]
        vtk["Momentum"] = @views (plot_gridmomentum[:, 1], plot_gridmomentum[:, 2], zero_gp_array)
        vtk["Force"] = @views (plot_gridforce[:, 1], plot_gridforce[:, 2], zero_gp_array)
        pvd_grid[t] = vtk
    end

    return nothing
end
