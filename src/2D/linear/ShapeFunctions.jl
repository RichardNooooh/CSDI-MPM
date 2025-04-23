
function linearFunction(particle_pos::Vec2{Float64}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Float64
    distance::Vec2{Float64} = abs.(particle_pos - nodal_pos)
    if sum(distance .> nodal_width) > 0.0
        return 0.0
    end

    result::Vec2{Float64} = Vec2{Float64}(1.0, 1.0) - distance ./ nodal_width
    return result[1] * result[2]
end

function topHatConstantFunction(particle_pos::Vec2{Float64}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Float64
    distance::Vec2{Float64} = abs.(particle_pos - nodal_pos)
    if sum(distance .> nodal_width) > 0.0
        return 0.0
    end

    return 0.25
end

# Implements a linear shape function
# \phi(\mathbf{x})_I = \Big(1 - \dfrac{|x^1 - x_I^1|}{h_I^1}\Big)\Big(1 - \dfrac{|x^2 - x_I^2|}{h_I^2}\Big)
# where x^j are the components of \mathbf{x}
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

function linearGradient(particle_pos::Vec2{Float64}, nodal_pos::Vec2{Float64}, nodal_width::Vec2{Float64})::Vec2{Float64}
    distance::Vec2{Float64} = particle_pos - nodal_pos
    abs_distance::Vec2{Float64} = abs.(distance)
    if sum(abs_distance .> nodal_width) > 0.0
        return 0.0
    end

    denominator::Float64 = nodal_width[1] * nodal_width[2]
    grad_1::Float64 = sign(distance[1]) * (abs_distance[2] - nodal_width[2]) / denominator
    grad_2::Float64 = sign(distance[2]) * (abs_distance[1] - nodal_width[1]) / denominator
    gradient::Vec2{Float64} = Vec2{Float64}(grad_1, grad_2)

    return gradient
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

function getQ4Centroid(material_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
    sum_vertices::Vec2{Float64} = sum(material_vertices)
    return sum_vertices / 4
end

function getL2Centroid(material_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
    sum_vertices::Vec2{Float64} = sum(material_vertices)
    return sum_vertices / 2
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
function findConnectedGrid(position::Vec2{Float64}, grid_cell_length::Vec2{Float64}, grid_num_nodes::Vec2{Int64})::Vec4{UInt32}
    botleft_grid::Vec2{Int64} = Int64.(floor.(position ./ grid_cell_length)) + Vec2{Int64}(1, 1)

    @assert botleft_grid[1] >= 1 && botleft_grid[1] < grid_num_nodes[1] string("Vertex found outside of the grid with index: ", botleft_grid[1], " from ", position[1])
    @assert botleft_grid[2] >= 1 && botleft_grid[2] < grid_num_nodes[2] string("Vertex found outside of the grid with index: ", botleft_grid[2], " from ", position[2])
    
    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]

    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    
    return Vec4{UInt32}(Index2DTo1D(i    , j    , numHorizontalGridNodes),
                        Index2DTo1D(i    , j + 1, numHorizontalGridNodes),
                        Index2DTo1D(i + 1, j    , numHorizontalGridNodes),
                        Index2DTo1D(i + 1, j + 1, numHorizontalGridNodes))
end

function getCellIndex(position::Vec2{Float64}, grid_num_nodes::Vec2{Int64}, grid_cell_length::Vec2{Float64})::UInt32
    botleft_grid::Vec2{Int64} = Int64.(floor.(position ./ grid_cell_length)) + Vec2{Int64}(1, 1)
    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]
    return Index2DTo1D(i    , j    , numHorizontalGridNodes)
end


