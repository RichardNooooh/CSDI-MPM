
const Point = Union{MaterialPoint, SurfacePoint, RigidPoint}

"""
Fills in the lookup table with coordinates. Note that vertices will move each time step
"""
function fillVertexCoordinates!(lookupTable::VTKLookup, domain::Vector{<:Point})::Nothing
    @inbounds for p in domain
        for (vertex_idx, vtk_tag) in enumerate(p.node_connectivity_vtk)
            pos = lookupTable.to_pos[vtk_tag]
            p_vertex = p.vertices[vertex_idx]
            lookupTable.x[pos] = p_vertex[1]
            lookupTable.y[pos] = p_vertex[2]
            lookupTable.z[pos] = p_vertex[3]
        end
    end
    return nothing
end

"""
Constructs computational domains and VTK lookup tables with input GMSH file
"""
function initialize3DGmshBulkAndSurface!(file_name::String)::Tuple{Vector{MaterialPoint}, Vector{SurfacePoint}, VTKLookup, VTKLookup}
    gmsh.initialize()
    gmsh.open(file_name)

    # Get global node table
    all_tags, all_vertices, _ = gmsh.model.mesh.getNodes(-1, -1, true, false)
    coords = reinterpret(Vec3{Float64}, all_vertices)
    tag2idx = Dict(all_tags .=> eachindex(all_tags))

    material_domain::Vector{MaterialPoint}, bulk_lookup::VTKLookup = parseBulk!(coords, tag2idx)
    surface_domain::Vector{SurfacePoint}, surf_lookup::VTKLookup = parseSurface!(coords, tag2idx)

    gmsh.finalize()

    return material_domain, surface_domain, bulk_lookup, surf_lookup
end

function parseBulk!(coords::AbstractArray{Vec3{Float64}}, tag2idx::Dict{UInt64, Int64})::Tuple{Vector{MaterialPoint}, VTKLookup}
    # material points are defined at each quadrangle element on the "Bulk" physical group (see msh file)
    # assume that this physical group is on physical tag 1 (dimension 2)
    tet_type::Int64 = gmsh.model.mesh.getElementType("Tetrahedron", 1)
    bulk_entities::Vector{Int64} = gmsh.model.getEntitiesForPhysicalGroup(3, 1)

    # first pass to construct bulk lookup table
    bulk_node_tags::Vector{Int64} = Vector{Int64}(undef, 0)
    material_domain::Vector{MaterialPoint} = Vector{MaterialPoint}(undef, 0)

    for entity in bulk_entities
        element_node_tags::Vector{Int64}, _, _ = gmsh.model.mesh.getNodesByElementType(tet_type, entity)
        append!(bulk_node_tags, element_node_tags)

        n_tets::Int64 = length(element_node_tags) ÷ 4
        sizehint!(material_domain, length(material_domain) + n_tets)

        for i in 1:4:length(element_node_tags)
            tag1, tag2, tag3, tag4 = element_node_tags[i:i+3]

            material_point::MaterialPoint = MaterialPoint()
            material_point.node_connectivity_vtk = Vec4{Int64}(tag1, tag2, tag3, tag4)
            push!(material_domain, material_point)
        end
    end
    bulk_lookup::VTKLookup = VTKLookup(unique(bulk_node_tags))

    # second pass to fill coordinates and geometry quantities
    @threads for material_point::MaterialPoint in material_domain
        vertices::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, 4)
        for (vertex_idx, vtk_tag) in enumerate(material_point.node_connectivity_vtk)
            vtk_vertex::Vec3{Float64} = coords[tag2idx[vtk_tag]]
            vertices[vertex_idx] = Vec3{Float64}(vtk_vertex[1], vtk_vertex[2], vtk_vertex[3])
        end
        material_point.vertices = vertices
        material_point.volume = material_point.volume_init = getTet4Volume(vertices)
    end

    return material_domain, bulk_lookup
end

function parseSurface!(coords::AbstractArray{Vec3{Float64}}, tag2idx::Dict{UInt64, Int64})::Tuple{Vector{SurfacePoint}, VTKLookup}
    tri_type::Int64 = gmsh.model.mesh.getElementType("Triangle", 1)
    surface_entities::Vector{Int64} = gmsh.model.getEntitiesForPhysicalGroup(2, 2)

    # first pass to construct bulk lookup table
    surface_node_tags::Vector{Int64} = Vector{Int64}(undef, 0)
    surface_domain::Vector{SurfacePoint} = Vector{SurfacePoint}(undef, 0)

    for entity in surface_entities
        element_node_tags::Vector{Int64}, _, _ = gmsh.model.mesh.getNodesByElementType(tri_type, entity)
        append!(surface_node_tags, element_node_tags)

        n_tris::Int64 = length(element_node_tags) ÷ 3
        sizehint!(surface_domain, length(surface_domain) + n_tris)

        for i in 1:3:length(element_node_tags)
            tag1, tag2, tag3 = element_node_tags[i:i+2]

            surface_point::SurfacePoint = SurfacePoint()
            surface_point.node_connectivity_vtk = Vec3{Int64}(tag1, tag2, tag3)
            push!(surface_domain, surface_point)
        end
    end
    surface_lookup::VTKLookup = VTKLookup(unique(surface_node_tags))

    # second pass to fill in coordinates and geometric quantities
    @threads for surface_point::SurfacePoint in surface_domain
        vertices::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, 3)
        for (vertex_idx, vtk_tag) in enumerate(surface_point.node_connectivity_vtk)
            vtk_vertex::Vec3{Float64} = coords[tag2idx[vtk_tag]]
            vertices[vertex_idx] = Vec3{Float64}(vtk_vertex[1], vtk_vertex[2], vtk_vertex[3])
        end
        surface_point.vertices = vertices
        surface_point.area_init, surface_point.normal_init = getTri3AreaAndNormal(surface_point.vertices)
        surface_point.area, surface_point.normal = surface_point.area_init, surface_point.normal_init
    end

    return surface_domain, surface_lookup
end


function parseRigidBulk!(coords::AbstractArray{Vec3{Float64}}, tag2idx::Dict{UInt64, Int64})::Tuple{Vector{RigidPoint}, VTKLookup}
    # material points are defined at each quadrangle element on the "Bulk" physical group (see msh file)
    # assume that this physical group is on physical tag 1 (dimension 2)
    tet_type::Int64 = gmsh.model.mesh.getElementType("Tetrahedron", 1)
    bulk_entities::Vector{Int64} = gmsh.model.getEntitiesForPhysicalGroup(3, 1)

    # first pass to construct bulk lookup table
    bulk_node_tags::Vector{Int64} = Vector{Int64}(undef, 0)
    rigid_domain::Vector{RigidPoint} = Vector{RigidPoint}(undef, 0)

    for entity in bulk_entities
        element_node_tags::Vector{Int64}, _, _ = gmsh.model.mesh.getNodesByElementType(tet_type, entity)
        append!(bulk_node_tags, element_node_tags)

        n_tets::Int64 = length(element_node_tags) ÷ 4
        sizehint!(material_domain, length(material_domain) + n_tets)

        for i in 1:4:length(element_node_tags)
            tag1, tag2, tag3, tag4 = element_node_tags[i:i+3]

            rigid_point::RigidPoint = RigidPoint()
            rigid_point.node_connectivity_vtk = Vec4{Int64}(tag1, tag2, tag3, tag4)
            push!(rigid_domain, rigid_point)
        end
    end
    rigid_lookup::VTKLookup = VTKLookup(unique(bulk_node_tags))

    # second pass to fill coordinates and geometry quantities
    @threads for rigid_point::RigidPoint in rigid_domain
        vertices::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, 4)
        for (vertex_idx, vtk_tag) in enumerate(rigid_point.node_connectivity_vtk)
            vtk_vertex::Vec3{Float64} = coords[tag2idx[vtk_tag]]
            vertices[vertex_idx] = Vec3{Float64}(vtk_vertex[1], vtk_vertex[2], vtk_vertex[3])
        end
        rigid_point.vertices = vertices
    end

    return rigid_domain, rigid_lookup
end


"""
VTK Output data structures and writers
"""

struct BulkVTKData
    pressure::Array{Float64, 1}
    volume::Array{Float64, 1}
    volumeVerts::Array{Float64, 1}
    density::Array{Float64, 1}
    jacobian::Array{Float64, 1}
    velocity::Array{Float64, 2}
    stress::Array{Float64, 3}
    deformgrad::Array{Float64, 3}
end


struct SurfaceVTKData
    surface_tension::Array{Float64, 1}
    area::Array{Float64, 1}
    normal::Array{Float64, 2}
end


struct GridVTKData
    mass::Array{Float64, 1}
    momentum::Array{Float64, 2}
    force::Array{Float64, 2}
end

function initializeVTKBulkData(material_domain::Vector{MaterialPoint}, destination_folder::String)::Tuple{BulkVTKData, WriteVTK.CollectionFile, Vector{MeshCell}}
    numMP::Int64 = length(material_domain)
    bulk_data = BulkVTKData(
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP, 3),
        Array{Float64}(undef, numMP, 3, 3),
        Array{Float64}(undef, numMP, 3, 3)
    )
    pvd_bulk = paraview_collection(joinpath(destination_folder, "bulk"))
    vtk_bulk_cells = [MeshCell(VTKCellTypes.VTK_TETRA, material_domain[mp_index].node_connectivity_vtk) for mp_index = 1:numMP]

    @info("Allocated $(Base.summarysize(bulk_data) + Base.summarysize(vtk_bulk_cells)) bytes for plotting the bulk")

    return bulk_data, pvd_bulk, vtk_bulk_cells
end

function initializeVTKSurfaceData(surface_domain::Vector{SurfacePoint}, destination_folder::String)::Tuple{SurfaceVTKData, WriteVTK.CollectionFile, Vector{MeshCell}}
    numSP::Int64 = length(surface_domain)
    surface_data = SurfaceVTKData(
        Array{Float64}(undef, numSP),
        Array{Float64}(undef, numSP),
        Array{Float64}(undef, numSP, 3)
    )
    pvd_surf = paraview_collection(joinpath(destination_folder, "surf"))
    vtk_surf_cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, surface_domain[sp_index].node_connectivity_vtk) for sp_index = 1:numSP]

    @info("Allocated $(Base.summarysize(surface_data) + Base.summarysize(vtk_surf_cells)) bytes for plotting the surface")

    return surface_data, pvd_surf, vtk_surf_cells
end


function initializeVTKGridData(grid::Grid, destination_folder::String)::Tuple{GridVTKData, WriteVTK.CollectionFile}
    numGP::Int64 = length(grid.points)
    grid_data = GridVTKData(
        Array{Float64}(undef, numGP),
        Array{Float64}(undef, numGP, 3),
        Array{Float64}(undef, numGP, 3)
    )
    pvd_grid = paraview_collection(joinpath(destination_folder, "grid"))

    @info("Allocated $(Base.summarysize(grid_data)) bytes for plotting the grid")

    return grid_data, pvd_grid
end


function initializeVTKRigidData(rigid_domain::Vector{RigidPoint}, destination_folder::String)::Tuple{WriteVTK.CollectionFile, Vector{MeshCell}}
    numRP::Int64 = length(rigid_domain)
    pvd_rigid = paraview_collection(joinpath(destination_folder, "rigid"))
    vtk_rigid_cells = [MeshCell(VTKCellTypes.VTK_TETRA, rigid_domain[rp_index].node_connectivity_vtk) for rp_index = 1:numRP]
    @info("Allocated $(Base.summarysize(vtk_rigid_cells)) bytes for plotting the rigid bulk")

    return pvd_rigid, vtk_rigid_cells
end


function writeVTIGrid(grid::Grid, destination_folder::String)::Nothing
    @info("Creating VTI file for the background grid")
    vtk_grid(joinpath(destination_folder, "grid"),
        0.0:grid.cell_length[1]:grid.grid_length[1],
        0.0:grid.cell_length[2]:grid.grid_length[2],
        0.0:grid.cell_length[3]:grid.grid_length[3]) do vtk
    end
    return nothing
end



function writeVTKBulk(material_domain::Vector{MaterialPoint},
                                bulk_data::BulkVTKData, bulk_vtk_lookup::VTKLookup,
                                vtk_bulk_cells::Vector{MeshCell}, pvd_bulk::WriteVTK.CollectionFile,
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    numMP::Int64 = length(material_domain)
    @threads for mp_index = 1:numMP
        material_point::MaterialPoint = material_domain[mp_index]

        # record all plotting information into these arrays
        bulk_data.pressure[mp_index] = material_point.pressure
        bulk_data.volume[mp_index] = material_point.volume
        bulk_data.volumeVerts[mp_index] = getTet4Volume(material_point.vertices)
        bulk_data.density[mp_index] = material_point.density
        bulk_data.jacobian[mp_index] = det(material_point.deformation_gradient)
        bulk_data.velocity[mp_index, 1] = material_point.velocity[1]
        bulk_data.velocity[mp_index, 2] = material_point.velocity[2]
        bulk_data.velocity[mp_index, 3] = material_point.velocity[3]

        bulk_data.stress[mp_index, 1, 1] = material_point.stress[1, 1]
        bulk_data.stress[mp_index, 1, 2] = material_point.stress[1, 2]
        bulk_data.stress[mp_index, 1, 3] = material_point.stress[1, 3]
        bulk_data.stress[mp_index, 2, 1] = material_point.stress[2, 1]
        bulk_data.stress[mp_index, 2, 2] = material_point.stress[2, 2]
        bulk_data.stress[mp_index, 2, 3] = material_point.stress[2, 3]
        bulk_data.stress[mp_index, 3, 1] = material_point.stress[3, 1]
        bulk_data.stress[mp_index, 3, 2] = material_point.stress[3, 2]
        bulk_data.stress[mp_index, 3, 3] = material_point.stress[3, 3]

        bulk_data.deformgrad[mp_index, 1, 1] = material_point.deformation_gradient[1, 1]
        bulk_data.deformgrad[mp_index, 1, 2] = material_point.deformation_gradient[1, 2]
        bulk_data.deformgrad[mp_index, 1, 3] = material_point.deformation_gradient[1, 3]
        bulk_data.deformgrad[mp_index, 2, 1] = material_point.deformation_gradient[2, 1]
        bulk_data.deformgrad[mp_index, 2, 2] = material_point.deformation_gradient[2, 2]
        bulk_data.deformgrad[mp_index, 2, 3] = material_point.deformation_gradient[2, 3]
        bulk_data.deformgrad[mp_index, 3, 1] = material_point.deformation_gradient[3, 1]
        bulk_data.deformgrad[mp_index, 3, 2] = material_point.deformation_gradient[3, 2]
        bulk_data.deformgrad[mp_index, 3, 3] = material_point.deformation_gradient[3, 3]
    end

    fillVertexCoordinates!(bulk_vtk_lookup, material_domain)

    # create vtk file
    vtk_grid(joinpath(vtk_output_folder, "bulk_timestep_$(Int64(t_index / t_record))"), 
                    bulk_vtk_lookup.x, bulk_vtk_lookup.y, bulk_vtk_lookup.z, vtk_bulk_cells) do vtk
        vtk["Pressure", VTKCellData()] =             bulk_data.pressure[:]
        vtk["Volume", VTKCellData()] =               bulk_data.volume[:]
        vtk["Volume From Vertices", VTKCellData()] = bulk_data.volumeVerts[:]
        vtk["Density", VTKCellData()] =              bulk_data.density[:]
        vtk["Jacobian", VTKCellData()] =             bulk_data.jacobian[:]
        vtk["Velocity", VTKCellData()] =             @views (bulk_data.velocity[:, 1], bulk_data.velocity[:, 2], bulk_data.velocity[:, 3])
        vtk["Stress", VTKCellData()] =               @views (bulk_data.stress[:, 1, 1], bulk_data.stress[:, 2, 1], bulk_data.stress[:, 3, 1],
                                                     bulk_data.stress[:, 1, 2], bulk_data.stress[:, 2, 2], bulk_data.stress[:, 3, 2],
                                                     bulk_data.stress[:, 1, 3], bulk_data.stress[:, 2, 3], bulk_data.stress[:, 3, 3])
        vtk["Deformation Gradient", VTKCellData()] = @views (bulk_data.deformgrad[:, 1, 1], bulk_data.deformgrad[:, 2, 1], bulk_data.deformgrad[:, 3, 1],
                                                     bulk_data.deformgrad[:, 1, 2], bulk_data.deformgrad[:, 2, 2], bulk_data.deformgrad[:, 3, 2],
                                                     bulk_data.deformgrad[:, 1, 3], bulk_data.deformgrad[:, 2, 3], bulk_data.deformgrad[:, 3, 3])
        pvd_bulk[t] = vtk
    end

    return nothing
end

function writeVTKSurface(surface_domain::Vector{SurfacePoint},
                                surface_data::SurfaceVTKData, surface_vtk_lookup::VTKLookup,
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

    vtk_grid(joinpath(vtk_output_folder, "surf_timestep_$(Int64(t_index / t_record))"), 
                    surface_vtk_lookup.x, surface_vtk_lookup.y, surface_vtk_lookup.z, vtk_surf_cells) do vtk
        vtk["Surface Tension", VTKCellData()]      = surface_data.surface_tension[:]
        vtk["Area", VTKCellData()]                 = surface_data.area[:]
        vtk["Normal", VTKCellData()]               = @views (surface_data.normal[:, 1], surface_data.normal[:, 2], surface_data.normal[:, 3])
        pvd_surf[t] = vtk
    end

    return nothing
end


function writeVTKGrid(grid::Grid,
                        grid_data::GridVTKData,
                        pvd_grid::WriteVTK.CollectionFile,
                        t_index::Int64, t_record::Int64, t::Float64,
                        vtk_output_folder::String)::Nothing
    numGP::Int64 = length(grid.points)
    for gp_index = 1:numGP
        grid_point::GridPoint = grid.points[gp_index]

        # record all plotting information into these arrays
        grid_data.mass[gp_index] = grid_point.mass
        grid_data.momentum[gp_index, 1] = grid_point.momentum[1]
        grid_data.momentum[gp_index, 2] = grid_point.momentum[2]
        grid_data.momentum[gp_index, 3] = grid_point.momentum[3]
        grid_data.force[gp_index, 1] = grid_point.force[1]
        grid_data.force[gp_index, 2] = grid_point.force[2]
        grid_data.force[gp_index, 3] = grid_point.force[3]
    end

    # create vtk file
    vtk_grid(joinpath(vtk_output_folder, "grid_timestep_$(Int64(t_index / t_record))"), 
             0.0:grid.cell_length[1]:grid.grid_length[1],
             0.0:grid.cell_length[2]:grid.grid_length[2]) do vtk
        vtk["Mass"] = grid_data.mass[:]
        vtk["Momentum"] = @views (grid_data.momentum[:, 1], grid_data.momentum[:, 2], grid_data.momentum[:, 3])
        vtk["Force"] = @views (grid_data.force[:, 1], grid_data.force[:, 2], grid_data.force[:, 3])
        pvd_grid[t] = vtk
    end

    return nothing
end


function writeVTKRigidBulk(rigid_domain::Vector{RigidPoint},
                                rigid_vtk_lookup::VTKLookup, vtk_rigid_cells::Vector{MeshCell},
                                pvd_rigid::WriteVTK.CollectionFile,
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    fillVertexCoordinates!(rigid_vtk_lookup, rigid_domain)

    vtk_grid(joinpath(vtk_output_folder, "rigid_timestep_$(Int64(t_index / t_record))"), 
                    rigid_vtk_lookup.x, rigid_vtk_lookup.y, rigid_vtk_lookup.z, vtk_rigid_cells) do vtk
        pvd_rigid[t] = vtk
    end
    return nothing
end
