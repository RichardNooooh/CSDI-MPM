
function initialize2DGmshBulkAndSurface!(file_name::String)::Vector{MaterialPoint}
    gmsh.initialize()
    gmsh.open(file_name)

    # for writing VTK files
    vtk_bulk_node_tags::Vector{Int64} = sort(unique(gmsh.model.mesh.getNodes(2, -1, true, false)[1])) # node tags in bulk
    # vtk_surf_node_tags::Vector{Int64} = sort(unique(gmsh.model.mesh.getNodes(1, -1, true, false)[1])) # node tags in surface

    # assume that this physical group is on physical tag 1 (dimension 2)
    material_domain::Vector{MaterialPoint} = parseBulk!(1, vtk_bulk_node_tags)

    # surface points are defined at each line element along the "Surface" physical group (see msh file)
    # assume that this physical group is on physical tag 2 (dimension 1)
    # surface_domain::Vector{SurfacePoint} = parseSurface!(2, vtk_surf_node_tags)

    # for writing VTKs
    gmsh.finalize()

    return material_domain
end

# function initialize2DGmshBilayer!(file_name::String, material_domain::Vector{MaterialPoint}, 
#                                                 monolayer_domain::Vector{SurfacePoint},
#                                                 bilayer_domain::Vector{SurfacePoint}, 
#                                                 solid_interface_domain::Vector{SurfacePoint},
#                                                 grid::Grid)::Nothing
#     gmsh.initialize()
#     gmsh.open(file_name)

#     # for writing VTK files
#     vtk_bulk_node_tags::Vector{Int64} = sort(unique(gmsh.model.mesh.getNodes(2, -1, true, false)[1])) # node tags in bulk
#     vtk_surf_node_tags::Vector{Int64} = sort(unique(gmsh.model.mesh.getNodes(1, -1, true, false)[1])) # node tags in surface

#     # assume that this physical group is on physical tag 1 (dimension 2)
#     parseBulk!(material_domain, grid, vtk_bulk_node_tags)

#     # surface points are defined at each line element along the "Surface" physical group (see msh file)
#     # assume that this physical group is on physical tag 2 through 4... (dimension 1)
#     parseSurface!(monolayer_domain, 2, grid, vtk_surf_node_tags)
#     parseSurface!(bilayer_domain, 3, grid, vtk_surf_node_tags)
#     parseSurface!(solid_interface_domain, 4, grid, vtk_surf_node_tags)

#     # for writing VTKs
#     gmsh.finalize()

#     return nothing
# end

function parseBulk!(physical_tag::Int64, vtk_bulk_node_tags::Array)::Vector{MaterialPoint}
    # material points are defined at each quadrangle element on the "Bulk" physical group (see msh file)
    # assume that this physical group is on physical tag 1 (dimension 2)
    quadType = gmsh.model.mesh.getElementType("Quadrangle", 1)
    bulk_entities = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag)

    mp_index::UInt64 = 1
    material_domain::Vector{MaterialPoint} = Vector{MaterialPoint}(undef, 0)
    
    for bulk_entity_tag in bulk_entities
        node_tags::Vector{Int64}, coords::Vector{Float64}, _ = gmsh.model.mesh.getNodesByElementType(quadType, bulk_entity_tag)
        append!(material_domain, Vector{MaterialPoint}(undef, div(length(node_tags), 4)))
        for i in 1:4:length(node_tags)
            # these should be in counter-clockwise order
            node_tag_1::Int64 = node_tags[i]
            node_tag_2::Int64 = node_tags[i + 1]
            node_tag_3::Int64 = node_tags[i + 2]
            node_tag_4::Int64 = node_tags[i + 3]

            # * make sure that the vtk_tags list contains the same node tags
            @assert node_tag_1 in vtk_bulk_node_tags && 
                node_tag_2 in vtk_bulk_node_tags && 
                node_tag_3 in vtk_bulk_node_tags && 
                node_tag_4 in vtk_bulk_node_tags

            position_1::Vec2{Float64} = Vec2{Float64}(coords[3 * i - 2] >= 0.0 ? coords[3 * i - 2] : 0.0, 
                                                        coords[3 * i - 1] >= 0.0 ? coords[3 * i - 1] : 0.0)
            position_2::Vec2{Float64} = Vec2{Float64}(coords[3 * i + 1] >= 0.0 ? coords[3 * i + 1] : 0.0, 
                                                        coords[3 * i + 2] >= 0.0 ? coords[3 * i + 2] : 0.0)
            position_3::Vec2{Float64} = Vec2{Float64}(coords[3 * i + 4] >= 0.0 ? coords[3 * i + 4] : 0.0, 
                                                        coords[3 * i + 5] >= 0.0 ? coords[3 * i + 5] : 0.0)
            position_4::Vec2{Float64} = Vec2{Float64}(coords[3 * i + 7] >= 0.0 ? coords[3 * i + 7] : 0.0, 
                                                        coords[3 * i + 8] >= 0.0 ? coords[3 * i + 8] : 0.0)

            # * make sure that coords and getNode coords are consistent
            # node_1_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_1)
            # node_2_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_2)
            # node_3_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_3)
            # node_4_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_4)
            # @assert node_1_coords[1] == position_1[1] && node_1_coords[2] == position_1[2]
            # @assert node_2_coords[1] == position_2[1] && node_2_coords[2] == position_2[2]
            # @assert node_3_coords[1] == position_3[1] && node_3_coords[2] == position_3[2]
            # @assert node_4_coords[1] == position_4[1] && node_4_coords[2] == position_4[2]
            
            # create material point
            material_point::MaterialPoint = MaterialPoint()
            # material_point.vertices = [position_1, position_2, position_3, position_4]
            vertices::Vector{Vec2{Float64}} = [position_1, position_2, position_3, position_4]
            material_point.position = getQ4Centroid(vertices)
            material_point.volume = getQ4Area(vertices)
            material_point.volume_init = material_point.volume

            # side lengths
            approx_side_length::Float64 = sqrt(material_point.volume)
            material_point.side_lengths = Vec2{Float64}(approx_side_length, approx_side_length)
            material_point.side_lengths_init = Vec2{Float64}(approx_side_length, approx_side_length)
            # # * make sure none of them is a Nothing object
            # for index::Int64 in material_point.node_connectivity_vtk
            #     @assert !isnothing(index)
            # end

            material_domain[mp_index] = material_point
            mp_index += 1
        end
    end

    return material_domain
end

# function parseSurface!(surface_physical_tag::Int64, vtk_surf_node_tags::Array)::Vector{SurfacePoint}
#     lineType = gmsh.model.mesh.getElementType("Line", 1)
#     surface_entities = gmsh.model.getEntitiesForPhysicalGroup(1, surface_physical_tag)
    
#     sp_index::UInt64 = 1
#     surface_domain::Vector{SurfacePoint} = Vector{SurfacePoint}(undef, 0)
    
#     for surface_entity_tag in surface_entities
#         node_tags::Vector{Int64}, coords::Vector{Float64}, _ = gmsh.model.mesh.getNodesByElementType(lineType, surface_entity_tag)
#         append!(surface_domain, Vector{SurfacePoint}(undef, div(length(node_tags), 2)))
#         for i in 1:2:length(node_tags)
#             node_tag_a::Int64 = node_tags[i + 1] # to orientate the normal correctly
#             node_tag_b::Int64 = node_tags[i]

#             # * make sure that the vtk_tags list contains the same node tags
#             @assert node_tag_a in vtk_surf_node_tags && 
#                 node_tag_b in vtk_surf_node_tags

#             position_a = Vec2{Float64}(coords[3 * i + 1] >= 0.0 ? coords[3 * i + 1] : 0.0, 
#                                         coords[3 * i + 2] >= 0.0 ? coords[3 * i + 2] : 0.0)  
#             position_b = Vec2{Float64}(coords[3 * i - 2] >= 0.0 ? coords[3 * i - 2] : 0.0, 
#                                         coords[3 * i - 1] >= 0.0 ? coords[3 * i - 1] : 0.0)
            
#             # * make sure that coords and getNode coords are consistent
#             # node_a_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_a)
#             # node_b_coords::Vector{Float64}, _, _, _ = gmsh.model.mesh.getNode(node_tag_b)
#             # @assert node_a_coords[1] == position_a[1] && node_a_coords[2] == position_a[2]
#             # @assert node_b_coords[1] == position_b[1] && node_b_coords[2] == position_b[2]

#             # create surface point
#             surface_point::SurfacePoint = SurfacePoint()
#             vertices::Vector{Vec2{Float64}} = [position_a, position_b]
#             surface_point.area, surface_point.normal = getL2AreaAndNormal(vertices)
#             surface_point.area_init, surface_point.normal_init = surface_point.area, surface_point.normal
#             surface_point.position = getL2Centroid(vertices)
            
#             surface_domain[sp_index] = surface_point
#             sp_index += 1
#         end
#     end

#     return surface_domain
# end
