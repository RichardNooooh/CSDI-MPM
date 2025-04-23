
mutable struct MaterialPoint
    # Properties
    velocity::Vec3{Float64}
    force_external::Vec3{Float64}
    stress::Mat33{Float64}
    deformation_gradient::Mat33{Float64}
    relative_deformation_gradient::Mat33{Float64}
    velocity_gradient::Mat33{Float64}

    mass::Float64
    volume::Float64
    volume_init::Float64
    density::Float64
    pressure::Float64

    parameter_1::Float64
    parameter_2::Float64

    # VTK file writing
    node_connectivity_vtk::Vec4{Int64}

    # connectivity
    connected_grid_array_length::UInt8     # denotes the final valid index in `connected_grid_indices`
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 32

    # CPDI-Tet4 vertex positions
    vertices::Vector{Vec3{Float64}}        # one-time allocated Vector of length 4

    function MaterialPoint()
        new(ZERO_VEC3,                   # velocity
            ZERO_VEC3,                   # force_external
            ZERO_MAT33,                  # stress
            IDENTITY_MAT33,              # deformation_gradient
            IDENTITY_MAT33,              # relative_deformation_gradient
            ZERO_MAT33,                  # velocity_gradient
            0.0,                         # mass
            0.0,                         # volume
            0.0,                         # volume_init
            0.0,                         # density
            0.0,                         # pressure
            0.0,                         # parameter_1
            0.0,                         # parameter_2
            Vec4{Int64}(0, 0, 0, 0),     # node_connectivity_vtk
            0,                           # connected_grid_array_length
            zeros(UInt32, 32),           # connected_grid_indices
            [ZERO_VEC3,                  # vertices
             ZERO_VEC3,
             ZERO_VEC3,
             ZERO_VEC3])
    end
end
