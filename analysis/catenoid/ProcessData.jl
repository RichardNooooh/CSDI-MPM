using Roots

function calculatedExpectedCatenoidConstant()::Float64
    D::Float64 = 1.0
    R::Float64 = 2.0

    f(A) = A * cosh(D / A) - R
    expected_A = find_zero(f, (1.6, 1.8))
    return expected_A
end


function isActiveGridPoint(grid_point::GridPoint, arguments::Dict{String,Any})::Bool
    return sqrt(grid_point.position[1]^2 + grid_point.position[2]^2) >= arguments["grid_r_min_active"] &&
           sqrt(grid_point.position[1]^2 + grid_point.position[2]^2) <= arguments["grid_r_max_active"]
end

# Locates all unique vertices and return list of points in coordinates (z, r)
function recordRadiiPoints(membrane_domain::Vector{MembranePoint})::Vector{Vec3{Float64}}
    radius_list::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, length(membrane_domain))
    for sp_index in eachindex(membrane_domain)
        membrane_point::MembranePoint = membrane_domain[sp_index]

        centroid::Vec3{Float64} = sum(membrane_point.vertices) / 3.0
        z::Float64 = centroid[3]
        r::Float64 = sqrt(centroid[1]^2 + centroid[2]^2)
        a::Float64 = membrane_point.area

        radius_list[sp_index] = Vec3{Float64}(z, r, a)
    end

    return radius_list
end

function record_current_radial_positions(membrane_domain::Vector{MembranePoint}, initial_centroid_position::Array{Float64,2})::Nothing
    @threads for sp_index in eachindex(membrane_domain)
        material_point::MembranePoint = membrane_domain[sp_index]
        centroid::Vec3{Float64} = ZERO_VEC3
        for vertex::Vec3{Float64} in material_point.vertices
            centroid += vertex
        end
        centroid /= 3

        initial_centroid_position[sp_index, 1] = centroid[1]
        initial_centroid_position[sp_index, 2] = centroid[2]
        initial_centroid_position[sp_index, 3] = 0.0
    end
end

function catenoid(z::Float64, A::Float64)
    return A * cosh((z - 1.0) / A)
end

function writeVTKMembraneCustom(membrane_domain::Vector{MembranePoint},
                                membrane_data::MembraneVTKData, radii_diff::Array{Float64, 1}, expected_catenoid_A::Float64,
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
        r = hypot(membrane_vtk_lookup.x[i], membrane_vtk_lookup.y[i])
        expected_radius::Float64 = expected_catenoid_A * cosh((membrane_vtk_lookup.z[i] - 1.0) / expected_catenoid_A)
        radii_diff[i] = abs(r - expected_radius)
    end

    vtk_grid(joinpath(vtk_output_folder, "membrane_timestep_$(Int64(t_index / t_record))"), 
                    membrane_vtk_lookup.x, membrane_vtk_lookup.y, membrane_vtk_lookup.z, vtk_membrane_cells) do vtk
        vtk["Surface Tension", VTKCellData()] = membrane_data.surface_tension[:]
        vtk["Area", VTKCellData()]            = membrane_data.area[:]
        vtk["Normal", VTKCellData()]          = @views (membrane_data.normal[:, 1], membrane_data.normal[:, 2], membrane_data.normal[:, 3])
        vtk["Velocity", VTKCellData()]        = @views (membrane_data.velocity[:, 1], membrane_data.velocity[:, 2], membrane_data.velocity[:, 3])
        vtk["Absolute Radii Difference", VTKPointData()] = radii_diff[:]
        pvd_membrane[t] = vtk
    end

    return nothing
end

