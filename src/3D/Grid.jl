
mutable struct GridPoint
    # nodal properties
    position::Vec3{Float64}
    momentum::Vec3{Float64}
    momentum_next::Vec3{Float64}
    velocity::Vec3{Float64}
    velocity_next::Vec3{Float64}
    force::Vec3{Float64}
    mass::Float64

    # Dirichlet boundary condition
    is_fixed::Vec3{Bool}

    # locking for parallelizing p2g methods
    lock::SpinLock

    function GridPoint(pos::Vec3{Float64})
        new(pos,                                # position
            ZERO_VEC3,                          # momentum
            ZERO_VEC3,                          # momentum_next
            ZERO_VEC3,                          # velocity
            ZERO_VEC3,                          # velocity_next
            ZERO_VEC3,                          # force
            0.0,                                # mass
            Vec3{Bool}(false, false, false),    # is_fixed
            SpinLock())                         # lock
    end
end

function Index3DTo1D(i::Int64, j::Int64, k::Int64, num_rows::Int64, num_cols::Int64)::UInt32
    return UInt32(i + (j - 1) * num_rows + (k - 1) * num_rows * num_cols)
end

struct Grid
    # grid properties
    grid_length::Vec3{Float64}
    num_nodes::Vec3{Int64}
    cell_length::Vec3{Float64}

    points::Vector{GridPoint}

    function Grid(l_x::Float64, l_y::Float64, l_z::Float64, n_x::Int64, n_y::Int64, n_z::Int64)
        distance_x::Float64 = l_x / Float64(n_x - 1.0)
        distance_y::Float64 = l_y / Float64(n_y - 1.0)
        distance_z::Float64 = l_z / Float64(n_z - 1.0)
        points_vec::Vector{GridPoint} = Vector{GridPoint}(undef, n_x * n_y * n_z)

        for k = 1:n_z
            for j = 1:n_y
                for i = 1:n_x
                    pos_x::Float64 = (i - 1) * distance_x
                    pos_y::Float64 = (j - 1) * distance_y
                    pos_z::Float64 = (k - 1) * distance_z
                    
                    points_vec[Index3DTo1D(i, j, k, n_x, n_y)] = GridPoint(Vec3{Float64}(pos_x, pos_y, pos_z))
                end
            end
        end

        new(Vec3{Float64}(l_x, l_y, l_z),                       # grid_length
            Vec3{Int64}(n_x, n_y, n_z),                         # num_nodes
            Vec3{Float64}(distance_x, distance_y, distance_z),  # cell_length
            points_vec)                                         # points
    end
end

Base.show(io::IO, grid::Grid) = print(io, "Dimensions: $(grid.grid_length)\n",
                                          "Num Nodes: $(grid.num_nodes)\n",
                                          "Cell Size: $(grid.cell_length)")

function resetGrid!(grid::Grid)::Nothing
    @threads for grid_point::GridPoint in grid.points
        grid_point.mass          = 0.0
        grid_point.momentum      = ZERO_VEC3
        grid_point.force         = ZERO_VEC3
        grid_point.velocity      = ZERO_VEC3
        grid_point.velocity_next = ZERO_VEC3
    end
    
    return nothing
end
