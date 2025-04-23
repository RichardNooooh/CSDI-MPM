
const Point = Union{MembranePoint, RigidPoint}

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
function initialize3DGmshSurface!(file_name::String)::Tuple{Vector{MembranePoint}, VTKLookup}
    gmsh.initialize()
    gmsh.open(file_name)

    # Get global node table
    all_tags, all_vertices, _ = gmsh.model.mesh.getNodes(-1, -1, true, false)
    coords = reinterpret(Vec3{Float64}, all_vertices)
    tag2idx = Dict(all_tags .=> eachindex(all_tags))

    membrane_domain::Vector{MembranePoint}, membrane_lookup::VTKLookup = parseSurface!(coords, tag2idx)

    gmsh.finalize()

    return membrane_domain, membrane_lookup
end


function parseSurface!(coords::AbstractArray{Vec3{Float64}}, tag2idx::Dict{UInt64, Int64})::Tuple{Vector{MembranePoint}, VTKLookup}
    tri_type::Int64 = gmsh.model.mesh.getElementType("Triangle", 1)
    surface_entities::Vector{Int64} = gmsh.model.getEntitiesForPhysicalGroup(2, 1)

    # first pass to construct bulk lookup table
    surface_node_tags::Vector{Int64} = Vector{Int64}(undef, 0)
    membrane_domain::Vector{MembranePoint} = Vector{MembranePoint}(undef, 0)

    for entity in surface_entities
        element_node_tags::Vector{Int64}, _, _ = gmsh.model.mesh.getNodesByElementType(tri_type, entity)
        append!(surface_node_tags, element_node_tags)

        n_tris::Int64 = length(element_node_tags) ÷ 3
        sizehint!(membrane_domain, length(membrane_domain) + n_tris)

        for i in 1:3:length(element_node_tags)
            tag1, tag2, tag3 = element_node_tags[i:i+2]

            membrane_point::MembranePoint = MembranePoint()
            membrane_point.node_connectivity_vtk = Vec3{Int64}(tag1, tag2, tag3)
            push!(membrane_domain, membrane_point)
        end
    end
    membrane_lookup::VTKLookup = VTKLookup(unique(surface_node_tags))

    # second pass to fill in coordinates and geometric quantities
    @threads for membrane_point::MembranePoint in membrane_domain
        vertices::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, 3)
        for (vertex_idx, vtk_tag) in enumerate(membrane_point.node_connectivity_vtk)
            vtk_vertex::Vec3{Float64} = coords[tag2idx[vtk_tag]]
            vertices[vertex_idx] = Vec3{Float64}(vtk_vertex[1], vtk_vertex[2], vtk_vertex[3])
        end
        membrane_point.vertices = vertices
        membrane_point.area, membrane_point.normal = getTri3AreaAndNormal(membrane_point.vertices)
    end

    return membrane_domain, membrane_lookup
end


"""
Constructs computational domains and VTK lookup tables with input GMSH file
"""
function initialize3DGmshRigidSurface!(file_name::String)::Tuple{Vector{RigidPoint}, VTKLookup}
    gmsh.initialize()
    gmsh.open(file_name)

    # Get global node table
    all_tags, all_vertices, _ = gmsh.model.mesh.getNodes(-1, -1, true, false)
    coords = reinterpret(Vec3{Float64}, all_vertices)
    tag2idx = Dict(all_tags .=> eachindex(all_tags))

    rigid_domain::Vector{RigidPoint}, rigid_lookup::VTKLookup = parseRigidSurface!(coords, tag2idx)

    gmsh.finalize()

    return rigid_domain, rigid_lookup
end


function parseRigidSurface!(coords::AbstractArray{Vec3{Float64}}, tag2idx::Dict{UInt64, Int64})::Tuple{Vector{RigidPoint}, VTKLookup}
    tri_type::Int64 = gmsh.model.mesh.getElementType("Triangle", 1)
    surface_entities::Vector{Int64} = gmsh.model.getEntitiesForPhysicalGroup(2, 1)

    # first pass to construct bulk lookup table
    rigid_node_tags::Vector{Int64} = Vector{Int64}(undef, 0)
    rigid_domain::Vector{RigidPoint} = Vector{RigidPoint}(undef, 0)

    for entity in surface_entities
        element_node_tags::Vector{Int64}, _, _ = gmsh.model.mesh.getNodesByElementType(tri_type, entity)
        append!(rigid_node_tags, element_node_tags)

        n_tris::Int64 = length(element_node_tags) ÷ 3
        sizehint!(rigid_domain, length(rigid_domain) + n_tris)

        for i in 1:3:length(element_node_tags)
            tag1, tag2, tag3 = element_node_tags[i:i+2]

            rigid_point::RigidPoint = RigidPoint()
            rigid_point.node_connectivity_vtk = Vec3{Int64}(tag1, tag2, tag3)
            push!(rigid_domain, rigid_point)
        end
    end
    rigid_lookup::VTKLookup = VTKLookup(unique(rigid_node_tags))

    # second pass to fill in coordinates and geometric quantities
    @threads for rigid_point::RigidPoint in rigid_domain
        vertices::Vector{Vec3{Float64}} = Vector{Vec3{Float64}}(undef, 3)
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

struct MembraneVTKData
    surface_tension::Array{Float64, 1}
    area::Array{Float64, 1}
    normal::Array{Float64, 2}
    velocity::Array{Float64, 2}
end

struct GridVTKData
    active::Array{Float64, 1}
    mass::Array{Float64, 1}
    momentum::Array{Float64, 2}
    force::Array{Float64, 2}
end


function initializeVTKMembraneData(membrane_domain::Vector{MembranePoint}, destination_folder::String)::Tuple{MembraneVTKData, WriteVTK.CollectionFile, Vector{MeshCell}}
    numMP::Int64 = length(membrane_domain)
    membrane_data = MembraneVTKData(
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP),
        Array{Float64}(undef, numMP, 3),
        Array{Float64}(undef, numMP, 3)
    )
    pvd_membrane = paraview_collection(joinpath(destination_folder, "membrane"))
    vtk_membrane_cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, membrane_domain[mp_index].node_connectivity_vtk) for mp_index = 1:numMP]

    @info("Allocated $(Base.summarysize(membrane_data) + Base.summarysize(vtk_membrane_cells)) bytes for plotting the membrane")

    return membrane_data, pvd_membrane, vtk_membrane_cells
end


function initializeVTKGridData(grid::Grid, destination_folder::String)::Tuple{GridVTKData, WriteVTK.CollectionFile}
    numGP::Int64 = length(grid.points)
    grid_data = GridVTKData(
        Array{Float64}(undef, numGP),
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
    vtk_rigid_cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, rigid_domain[rp_index].node_connectivity_vtk) for rp_index = 1:numRP]
    @info("Allocated $(Base.summarysize(vtk_rigid_cells)) bytes for plotting the rigid surface")

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


function writeVTKMembrane(membrane_domain::Vector{MembranePoint},
                                membrane_data::MembraneVTKData, membrane_vtk_lookup::VTKLookup,
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

    vtk_grid(joinpath(vtk_output_folder, "membrane_timestep_$(Int64(t_index / t_record))"), 
                    membrane_vtk_lookup.x, membrane_vtk_lookup.y, membrane_vtk_lookup.z, vtk_membrane_cells) do vtk
        vtk["Surface Tension", VTKCellData()] = membrane_data.surface_tension[:]
        vtk["Area", VTKCellData()]            = membrane_data.area[:]
        vtk["Normal", VTKCellData()]          = @views (membrane_data.normal[:, 1], membrane_data.normal[:, 2], membrane_data.normal[:, 3])
        vtk["Velocity", VTKCellData()]        = @views (membrane_data.velocity[:, 1], membrane_data.velocity[:, 2], membrane_data.velocity[:, 3])
        pvd_membrane[t] = vtk
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
             0.0:grid.cell_length[2]:grid.grid_length[2],
             0.0:grid.cell_length[3]:grid.grid_length[3]) do vtk
        vtk["Active"] = grid_data.active[:]
        vtk["Mass"] = grid_data.mass[:]
        vtk["Momentum"] = @views (grid_data.momentum[:, 1], grid_data.momentum[:, 2], grid_data.momentum[:, 3])
        vtk["Force"] = @views (grid_data.force[:, 1], grid_data.force[:, 2], grid_data.force[:, 3])
        pvd_grid[t] = vtk
    end

    return nothing
end


function writeVTKRigid(rigid_domain::Vector{RigidPoint}, rigid_vtk_lookup::VTKLookup,
                                vtk_rigid_cells::Vector{MeshCell}, pvd_rigid::WriteVTK.CollectionFile, 
                                t_index::Int64, t_record::Int64, t::Float64,
                                vtk_output_folder::String)::Nothing
    fillVertexCoordinates!(rigid_vtk_lookup, rigid_domain)

    vtk_grid(joinpath(vtk_output_folder, "rigid_timestep_$(Int64(t_index / t_record))"), 
                    rigid_vtk_lookup.x, rigid_vtk_lookup.y, rigid_vtk_lookup.z, vtk_rigid_cells) do vtk
        pvd_rigid[t] = vtk
    end

    return nothing
end
