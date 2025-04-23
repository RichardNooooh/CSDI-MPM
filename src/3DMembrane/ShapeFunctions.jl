
# Implements a linear shape function
# \phi(\mathbf{x})_I = \Big(1 - \dfrac{|x^1 - x_I^1|}{h_I^1}\Big)\Big(1 - \dfrac{|x^2 - x_I^2|}{h_I^2}\Big)\Big(1 - \dfrac{|x^3 - x_I^3|}{h_I^3}\Big)
# where x^j are the components of \mathbf{x}
function linearFunction(particle_pos::Vec3{Float64}, nodal_pos::Vec3{Float64}, nodal_width::Vec3{Float64})::Float64
    distance::Vec3{Float64} = abs.(particle_pos - nodal_pos)
    if sum(distance .> nodal_width) > 0.0
        return 0.0
    end

    result::Vec3{Float64} = Vec3{Float64}(1.0, 1.0, 1.0) - distance ./ nodal_width
    return result[1] * result[2] * result[3]
end

function surfCPDITri3ShapeFunctionAndGradient(surface_vertices::Vector{Vec3{Float64}}, nodal_pos::Vec3{Float64}, nodal_width::Vec3{Float64})::Tuple{Float64, Vec3{Float64}}
    # points
    p1::Vec3{Float64} = surface_vertices[1]; p2::Vec3{Float64} = surface_vertices[2]; p3::Vec3{Float64} = surface_vertices[3]

    # N_I
    linear1::Float64 = linearFunction(p1, nodal_pos, nodal_width)
    linear2::Float64 = linearFunction(p2, nodal_pos, nodal_width)
    linear3::Float64 = linearFunction(p3, nodal_pos, nodal_width)

    # geometric area (x2) and normal
    n::Vec3{Float64} = cross(p2 - p1, p3 - p1)
    J::Float64 = norm(n)
    n /= J
    scale_factor::Float64 = 1/J


    pc::Vec3{Float64} = (p1 + p2 + p3) / 3
    S_1::Float64 = dot(n, cross(p2, p3) + cross(pc, p2-p3))
    S_2::Float64 = dot(n, cross(p3, p1) + cross(pc, p3-p1))
    S_3::Float64 = dot(n, cross(p1, p2) + cross(pc, p1-p2))
    value::Float64 = scale_factor * (S_1*linear1 + S_2*linear2 + S_3*linear3)
    grad::Vec3{Float64} = scale_factor * cross((linear1*(p2-p3)) + (linear2*(p3-p1)) + (linear3*(p1-p2)), n)

    return value, grad
end

function getTri3Normal(surface_vertices::Vector{Vec3{Float64}})::Vec3{Float64}
    p1::Vec3{Float64} = surface_vertices[1]; p2::Vec3{Float64} = surface_vertices[2]; p3::Vec3{Float64} = surface_vertices[3]

    n::Vec3{Float64} = cross(p2 - p1, p3 - p1)
    n /= norm(n)
    return n
end

function getTri3Area(surface_vertices::Vector{Vec3{Float64}})::Float64
    p1::Vec3{Float64} = surface_vertices[1]; p2::Vec3{Float64} = surface_vertices[2]; p3::Vec3{Float64} = surface_vertices[3]

    n::Vec3{Float64} = cross(p2 - p1, p3 - p1)
    return norm(n) / 2
end

function getTri3AreaAndNormal(surface_vertices::Vector{Vec3{Float64}})::Tuple{Float64, Vec3{Float64}}
    # points
    p1::Vec3{Float64} = surface_vertices[1]; p2::Vec3{Float64} = surface_vertices[2]; p3::Vec3{Float64} = surface_vertices[3]
    
    n::Vec3{Float64} = cross(p2 - p1, p3 - p1)
    A::Float64 = norm(n) / 2
    n /= 2 * A

    return (A, n)
end

# finds the grid indices within the cell of a vertex_position
function findConnectedGrid(vertex_position::Vec3{Float64}, grid_cell_length::Vec3{Float64}, grid_num_nodes::Vec3{Int64})::Vec8{UInt32}
    botleft_grid::Vec3{Int64} = Int64.(floor.(vertex_position ./ grid_cell_length)) + Vec3{Int64}(1, 1, 1)

    @assert botleft_grid[1] >= 1 string("We got a negative grid[1] index: ", botleft_grid[1], " from ", vertex_position[1])
    @assert botleft_grid[2] >= 1 string("We got a negative grid[2] index: ", botleft_grid[2], " from ", vertex_position[2])
    @assert botleft_grid[3] >= 1 string("We got a negative grid[3] index: ", botleft_grid[3], " from ", vertex_position[3])

    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]
    k::Int64 = botleft_grid[3]

    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    numVerticalGridNodes::Int64 = grid_num_nodes[2]
    
    return Vec8{UInt32}(Index3DTo1D(i    , j    , k    , numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i    , j    , k + 1, numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i    , j + 1, k    , numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i    , j + 1, k + 1, numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i + 1, j    , k    , numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i + 1, j    , k + 1, numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i + 1, j + 1, k    , numHorizontalGridNodes, numVerticalGridNodes),
                        Index3DTo1D(i + 1, j + 1, k + 1, numHorizontalGridNodes, numVerticalGridNodes))
end


function getAllConnectedGrid!(connected_grid_indices::Vector{UInt32}, vertices::Vector{Vec3{Float64}}, grid_cell_length::Vec3{Float64}, grid_num_nodes::Vec3{Int64})::UInt8
    for vertex_index in 1:length(vertices)
        vertex::Vec3{Float64} = vertices[vertex_index]
        adjacent_grid_nodes::Vec8{UInt32} = findConnectedGrid(vertex, grid_cell_length, grid_num_nodes)

        maxConnectedNodes = 8
        for i = 1:maxConnectedNodes
            connected_grid_indices[maxConnectedNodes * vertex_index - (maxConnectedNodes - i)] = adjacent_grid_nodes[i]
        end
    end

    sort!(connected_grid_indices)
    return arrangeUniqueInPlace!(connected_grid_indices)
end
