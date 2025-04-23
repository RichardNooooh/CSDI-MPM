
function Index2DTo1D(i::Int64, j::Int64, num_rows::Int64)::UInt32
    return UInt32(i + (j - 1) * num_rows)
end

mutable struct GridPoint
    # nodal properties
    position::Vec2{Float64}
    momentum::Vec2{Float64}
    momentum_next::Vec2{Float64}
    velocity::Vec2{Float64}
    velocity_next::Vec2{Float64}
    force::Vec2{Float64}
    mass::Float64

    # Dirichlet boundary condition
    is_fixed::Vec2{Bool}

    # locking for parallelizing p2g methods
    lock::SpinLock

    function GridPoint(pos::Vec2{Float64})
        new(pos,                        # position
            ZERO_VEC2,                  # momentum
            ZERO_VEC2,                  # momentum_next
            ZERO_VEC2,                  # velocity
            ZERO_VEC2,                  # velocity_next
            ZERO_VEC2,    # force
            0.0,                        # mass
            Vec2{Bool}(false, false),   # is_fixed
            SpinLock())                 # lock
    end
end

struct Grid
    # grid properties
    grid_length::Vec2{Float64}
    num_nodes::Vec2{Int64}
    cell_length::Vec2{Float64}

    points::Vector{GridPoint}

    function Grid(l_x::Float64, l_y::Float64, n_x::Int64, n_y::Int64)
        distance_x::Float64 = l_x / Float64(n_x - 1.0)
        distance_y::Float64 = l_y / Float64(n_y - 1.0)
        points_vec::Vector{GridPoint} = Vector{GridPoint}(undef, n_x * n_y)

        for j = 1:n_y
            for i = 1:n_x
                pos_x::Float64 = (i - 1) * distance_x
                pos_y::Float64 = (j - 1) * distance_y
                points_vec[Index2DTo1D(i, j, n_x)] = GridPoint(Vec2{Float64}(pos_x, pos_y))
            end
        end

        new(Vec2{Float64}(l_x, l_y),                # grid_length
            Vec2{Int64}(n_x, n_y),                  # num_nodes
            Vec2{Float64}(distance_x, distance_y),  # cell_length
            points_vec)                             # points
    end
end

Base.show(io::IO, grid::Grid) = print(io, "Dimensions: $(grid.grid_length)\n",
                                          "Num Nodes: $(grid.num_nodes)\n",
                                          "Cell Size: $(grid.cell_length)")

function resetGrid!(grid::Grid)::Nothing
    @threads for grid_point::GridPoint in grid.points
        grid_point.mass = 0.0
        grid_point.momentum = ZERO_VEC2
        grid_point.force = ZERO_VEC2
    end

    return nothing
end

