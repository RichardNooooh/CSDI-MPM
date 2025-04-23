
# Implements a linear shape function
# \phi(\mathbf{x})_I = \Big(1 - \dfrac{|x^1 - x_I^1|}{h_I^1}\Big)\Big(1 - \dfrac{|x^2 - x_I^2|}{h_I^2}\Big)
# where x^j are the components of \mathbf{x}
function linearFunction(particle_pos::Vec2{Float64}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Float64
    distance::Vec2{Float64} = abs.(particle_pos - nodal_pos)
    if sum(distance .> nodal_width) > 0.0
        return 0.0
    end

    result::Vec2{Float64} = Vec2{Float64}(1.0, 1.0) - distance ./ nodal_width
    return result[1] * result[2]
end

function linearFunctionAndGradient(particle_pos::Vec2{Float64}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Tuple{Float64, Vec2{Float64}}
    distance::Vec2{Float64} = particle_pos - nodal_pos
    abs_distance::Vec2{Float64} = abs.(distance)
    if sum(abs_distance .> nodal_width) > 0.0
        return 0.0, ZERO_VEC2
    end

    result::Vec2{Float64} = Vec2{Float64}(1.0, 1.0) - abs_distance ./ nodal_width
    value::Float64 = result[1] * result[2]

    denominator::Float64 = nodal_width[1] * nodal_width[2]
    grad_1::Float64 = sign(distance[1]) * (abs_distance[2] - nodal_width[2]) / denominator
    grad_2::Float64 = sign(distance[2]) * (abs_distance[1] - nodal_width[1]) / denominator
    gradient::Vec2{Float64} = Vec2{Float64}(grad_1, grad_2)

    return value, gradient
end

function getQ4Area(material_vertices::Vector{Vec2{Float64}})::Float64
    p1::Vec2{Float64} = material_vertices[1]
    p2::Vec2{Float64} = material_vertices[2]
    p3::Vec2{Float64} = material_vertices[3]
    p4::Vec2{Float64} = material_vertices[4]

    pV::Float64 = 0.5 * ((p1[1] * p2[2] - p2[1] * p1[2])
                        + (p2[1] * p3[2] - p3[1] * p2[2])
                        + (p3[1] * p4[2] - p4[1] * p3[2])
                        + (p4[1] * p1[2] - p1[1] * p4[2]))

    return pV
end


# TODO documentation
# Page 53 of "Material Point Method after 25 Years" (Vaucorbeil, Nguyen, Sinaie, Wu)
function bulkCPDIShapeFunctionAndGradient(material_vertices::Vector{Vec2{Float64}}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Tuple{Float64, Vec2{Float64}}
    p1::Vec2{Float64} = material_vertices[1]
    p2::Vec2{Float64} = material_vertices[2]
    p3::Vec2{Float64} = material_vertices[3]
    p4::Vec2{Float64} = material_vertices[4]

    pV::Float64 = getQ4Area(material_vertices)

    linear1::Float64 = linearFunction(p1, nodal_pos, nodal_width)
    linear2::Float64 = linearFunction(p2, nodal_pos, nodal_width)
    linear3::Float64 = linearFunction(p3, nodal_pos, nodal_width)
    linear4::Float64 = linearFunction(p4, nodal_pos, nodal_width)

    # phi computation
    a::Float64 = (p4[1] - p1[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p4[2] - p1[2])
    b::Float64 = (p3[1] - p4[1]) * (p1[2] - p2[2]) - (p1[1] - p2[1]) * (p3[2] - p4[2])

    value1::Float64 = (6 * pV - a - b) * linear1
    value2::Float64 = (6 * pV - a + b) * linear2
    value3::Float64 = (6 * pV + a + b) * linear3
    value4::Float64 = (6 * pV + a - b) * linear4

    value::Float64 = (1 / (24 * pV)) * (value1 + value2 + value3 + value4)
    
    # gradient of phi computation
    grad1::Vec2{Float64} = Vec2{Float64}(p2[2] - p4[2], p4[1] - p2[1]) * linear1
    grad2::Vec2{Float64} = Vec2{Float64}(p3[2] - p1[2], p1[1] - p3[1]) * linear2
    grad3::Vec2{Float64} = Vec2{Float64}(p4[2] - p2[2], p2[1] - p4[1]) * linear3
    grad4::Vec2{Float64} = Vec2{Float64}(p1[2] - p3[2], p3[1] - p1[1]) * linear4

    gradient::Vec2{Float64} = (1 / (2 * pV)) * (grad1 + grad2 + grad3 + grad4)

    return value, gradient
end

function surfCPDIL2ShapeFunctionAndGradient(surface_vertices::Vector{Vec2{Float64}}, 
                            nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Tuple{Float64, Vec2{Float64}}
    x_a::Vec2{Float64} = surface_vertices[1]
    x_b::Vec2{Float64} = surface_vertices[2]
    
    A::Float64 = getL2Area(surface_vertices)
    scale::Float64 = 1 / A^2
    linear_a::Float64 = linearFunction(x_a, nodal_pos, nodal_width)
    linear_b::Float64 = linearFunction(x_b, nodal_pos, nodal_width)

    value::Float64 = 0.5*(linear_a + linear_b)
    gradient::Vec2{Float64} = scale * (linear_b - linear_a) * (x_b - x_a)
    return value, gradient
end

function getL2Area(surface_vertices::Vector{Vec2{Float64}})::Float64
    p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

    n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
    a::Float64 = norm(n)

    return a
end

function getL2Normal(surface_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
    p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

    n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
    n /= norm(n)
    return n
end

function getL2AreaAndNormal(surface_vertices::Vector{Vec2{Float64}})::Tuple{Float64, Vec2{Float64}}
    p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

    n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
    a::Float64 = norm(n)
    n /= a

    return a, n
end

# finds the grid indices within the cell of a vertex_position
function findConnectedGrid(vertex_position::Vec2{Float64}, grid_cell_length::Vec2{Float64}, grid_num_nodes::Vec2{Int64})::Vec4{UInt32}
    botleft_grid::Vec2{Int64} = Int64.(floor.(vertex_position ./ grid_cell_length)) + Vec2{Int64}(1, 1)

    @assert botleft_grid[1] >= 1 && botleft_grid[1] < grid_num_nodes[1] string("Vertex found outside of the grid with index: ", botleft_grid[1], " from ", vertex_position[1])
    @assert botleft_grid[2] >= 1 && botleft_grid[2] < grid_num_nodes[2] string("Vertex found outside of the grid with index: ", botleft_grid[2], " from ", vertex_position[2])
    
    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]

    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    
    return Vec4{UInt32}(Index2DTo1D(i    , j    , numHorizontalGridNodes),
                        Index2DTo1D(i    , j + 1, numHorizontalGridNodes),
                        Index2DTo1D(i + 1, j    , numHorizontalGridNodes),
                        Index2DTo1D(i + 1, j + 1, numHorizontalGridNodes))
end

function getAllConnectedGrid!(connected_grid_indices::Vector{UInt32}, vertices::Vector{Vec2{Float64}}, grid_cell_length::Vec2{Float64}, grid_num_nodes::Vec2{Int64})::UInt8
    for vertex_index in 1:length(vertices)
        vertex::Vec2{Float64} = vertices[vertex_index]
        adjacent_grid_nodes::Vec4{UInt64} = findConnectedGrid(vertex, grid_cell_length, grid_num_nodes)

        connected_grid_indices[4 * vertex_index - 3] = adjacent_grid_nodes[1]
        connected_grid_indices[4 * vertex_index - 2] = adjacent_grid_nodes[2]
        connected_grid_indices[4 * vertex_index - 1] = adjacent_grid_nodes[3]
        connected_grid_indices[4 * vertex_index    ] = adjacent_grid_nodes[4]
    end

    sort!(connected_grid_indices)
    return arrangeUniqueInPlace!(connected_grid_indices)
end
