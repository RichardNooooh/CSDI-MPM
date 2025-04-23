
function GIMPShapeFunctionAndGradient(particle_pos::Vec2{Float64},
    nodal_pos::Vec2{Float64},
    particle_lengths::Vec2{Float64},
    nodal_width::Vec2{Float64})::Tuple{Float64, Vec2{Float64}}

    function GIMPValueComponent(x::Float64, h::Float64, l::Float64)::Float64
        half_l::Float64 = 0.5 * l
        abs_x::Float64 = abs(x)
        if abs_x < half_l
            return 1.0 - (4.0 * x^2 + l^2) / (4 * h * l)
        elseif half_l <= abs_x <= h - half_l
            return 1 - abs_x / h
        else
            return (h + half_l - abs_x)^2 / (2 * h * l)
        end

        # if -h - l < x <= -h + l
        #     return (h + l + x)^2 / (4 * h * l)
        # elseif -h + l < x <= -l
        #     return 1 + x / h
        # elseif -l < x <= l
        #     return 1 - (x^2 + l^2)/(2 * h * l)
        # elseif l < x <= h - l
        #     return 1 - x / h
        # elseif h - l < x <= h + l
        #     return (h + l - x)^2 / (4 * h * l)
        # else
        #     return 0.0
        # end
    end

    function GIMPGradientComponent(x::Float64, h::Float64, l::Float64)::Float64
        half_l::Float64 = 0.5 * l
        abs_x::Float64 = abs(x)
        if abs_x < half_l
            return -8 * x / (4 * h * l)
        elseif half_l <= abs_x <= h - half_l
            return -(1.0 / h) * sign(x)
        else
            return -sign(x) * (h + half_l - abs_x) / (h * l)
        end
        # if -h - l < x <= -h + l
        #     return (h + l + x) / (2 * h * l)
        # elseif -h + l < x <= -l
        #     return 1 / h
        # elseif -l < x <= l
        #     return -x / (h * l)
        # elseif l < x <= h - l
        #     return -1/ h
        # elseif h - l < x <= h + l
        #     return -(h + l - x) / (2 * h * l)
        # else
        #     return 0.0
        # end
    end

    distance::Vec2{Float64} = particle_pos - nodal_pos
    abs_distance::Vec2{Float64} = abs.(distance)
    if any(abs_distance .> (nodal_width + 0.5.*particle_lengths))
        return 0.0, ZERO_VEC2
    end

    value::Float64 = (GIMPValueComponent(distance[1], nodal_width[1], particle_lengths[1]) *
                      GIMPValueComponent(distance[2], nodal_width[2], particle_lengths[2]))
    gradient::Vec2{Float64} = Vec2{Float64}(GIMPGradientComponent(distance[1], nodal_width[1], particle_lengths[1]),
                                        GIMPGradientComponent(distance[2], nodal_width[2], particle_lengths[2]))
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

function getQ4Centroid(material_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
    sum_vertices::Vec2{Float64} = sum(material_vertices)
    return sum_vertices / 4
end

# function getL2Centroid(material_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
#     sum_vertices::Vec2{Float64} = sum(material_vertices)
#     return sum_vertices / 2
# end

# function getL2Area(surface_vertices::Vector{Vec2{Float64}})::Float64
#     p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

#     n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
#     a::Float64 = norm(n)

#     return a
# end

# function getL2Normal(surface_vertices::Vector{Vec2{Float64}})::Vec2{Float64}
#     p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

#     n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
#     n /= norm(n)
#     return n
# end

# function getL2AreaAndNormal(surface_vertices::Vector{Vec2{Float64}})::Tuple{Float64, Vec2{Float64}}
#     p1::Vec2{Float64} = surface_vertices[1]; p2::Vec2{Float64} = surface_vertices[2]

#     n::Vec2{Float64} = Vec2{Float64}(p2[2] - p1[2], p1[1] - p2[1])
#     a::Float64 = norm(n)
#     n /= a

#     return a, n
# end

# finds the grid indices within the cell of a vertex_position
# TODO fix and generalize - currently only valid for sufficiently small l_p
function findConnectedGrid(position::Vec2{Float64}, side_lengths::Vec2{Float64}, grid_cell_length::Vec2{Float64}, grid_num_nodes::Vec2{Int64})::Vector{UInt32}
    botleft_grid::Vec2{Int64} = Int64.(floor.(position ./ grid_cell_length)) + Vec2{Int64}(1, 1)

    @assert botleft_grid[1] >= 1 && botleft_grid[1] < grid_num_nodes[1] string("Vertex found outside of the grid with index: ", botleft_grid[1], " from ", position[1])
    @assert botleft_grid[2] >= 1 && botleft_grid[2] < grid_num_nodes[2] string("Vertex found outside of the grid with index: ", botleft_grid[2], " from ", position[2])
    
    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]

    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    return Vector{UInt32}([Index2DTo1D(i    , j    , numHorizontalGridNodes),
                         Index2DTo1D(i    , j + 1, numHorizontalGridNodes),
                         Index2DTo1D(i + 1, j    , numHorizontalGridNodes),
                         Index2DTo1D(i + 1, j + 1, numHorizontalGridNodes),

                         Index2DTo1D(i - 1, j - 1, numHorizontalGridNodes),
                         Index2DTo1D(i    , j - 1, numHorizontalGridNodes),
                         Index2DTo1D(i + 1, j - 1, numHorizontalGridNodes),
                         Index2DTo1D(i + 2, j - 1, numHorizontalGridNodes),

                         Index2DTo1D(i - 1, j    , numHorizontalGridNodes),
                         Index2DTo1D(i + 2, j    , numHorizontalGridNodes),
                         Index2DTo1D(i - 1, j + 1, numHorizontalGridNodes),
                         Index2DTo1D(i + 2, j + 1, numHorizontalGridNodes),

                         Index2DTo1D(i - 1, j + 2, numHorizontalGridNodes),
                         Index2DTo1D(i    , j + 2, numHorizontalGridNodes),
                         Index2DTo1D(i + 1, j + 2, numHorizontalGridNodes),
                         Index2DTo1D(i + 2, j + 2, numHorizontalGridNodes)])
end

function getCellIndex(position::Vec2{Float64}, grid_num_nodes::Vec2{Int64}, grid_cell_length::Vec2{Float64})::UInt32
    botleft_grid::Vec2{Int64} = Int64.(floor.(position ./ grid_cell_length)) + Vec2{Int64}(1, 1)
    numHorizontalGridNodes::Int64 = grid_num_nodes[1]
    i::Int64 = botleft_grid[1]
    j::Int64 = botleft_grid[2]
    return Index2DTo1D(i    , j    , numHorizontalGridNodes)
end


