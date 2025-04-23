mutable struct RigidPoint
    velocity::Vec3{Float64}
    
    # Dirichlet Boundary Conditions
    is_fixed::Vec3{Bool}

    # VTK file writing
    node_connectivity_vtk::Vec3{Int64}

    # connectivity
    connected_grid_array_length::UInt8    # denotes the final valid index in `connected_grid_points`
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 24

    # CSDI-Tri3 vertex positions
    vertices::Vector{Vec3{Float64}}

    function RigidPoint()
        new(ZERO_VEC3,                    # velocity
            Vec3{Bool}(true, true, true), # is_fixed
            ZERO_VEC3,                    # node_connectivity_vtk
            0,                            # connected_grid_array_length
            zeros(UInt32, 24),            # connected_grid_points
            [ZERO_VEC3,                   # vertices
             ZERO_VEC3,
             ZERO_VEC3,
             ZERO_VEC3])
    end
end
