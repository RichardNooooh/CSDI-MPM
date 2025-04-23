mutable struct RigidPoint
    velocity::Vec2{Float64}

    # Dirichlet Boundary Conditions
    is_fixed::Vec2{Bool}

    # VTK file writing
    node_connectivity_vtk::Vec4{Int64}

    # connectivity
    connected_grid_array_length::UInt8    # denotes the final valid index in `connected_grid_points`
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 16

    # CPDI-Q4 vertex positions
    vertices::Vector{Vec2{Float64}}

    function RigidPoint()
        new(ZERO_VEC2,                 # velocity
            Vec2{Float64}(true, true), # is_fixed
            Vec4{Int64}(0, 0, 0, 0),   # node_connectivity_vtk
            0,                         # connected_grid_array_length
            zeros(UInt32, 16),         # connected_grid_points
            [ZERO_VEC2,                # vertices
             ZERO_VEC2,
             ZERO_VEC2,
             ZERO_VEC2])
    end
end
