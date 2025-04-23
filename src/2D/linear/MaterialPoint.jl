
mutable struct MaterialPoint
    # Properties
    position::Vec2{Float64}
    velocity::Vec2{Float64}
    force_external::Vec2{Float64}
    stress::Mat22{Float64}
    strain::Mat22{Float64}
    velocity_gradient::Mat22{Float64}

    deformation_gradient::Mat22{Float64}
    deformation_gradient_bar::Mat22{Float64}
    relative_deformation_gradient::Mat22{Float64}
    div_velocity::Float64

    mass::Float64
    volume::Float64
    volume_init::Float64
    density::Float64
    pressure::Float64

    parameter_1::Float64
    parameter_2::Float64

    bulk_ID::UInt8

    # connectivity
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 4

    function MaterialPoint()
        new(ZERO_VEC2,                        # position
            ZERO_VEC2,                        # velocity
            ZERO_VEC2,                        # force_external
            ZERO_MAT22,                                     # stress
            ZERO_MAT22,                                     # strain
            ZERO_MAT22,                                     # velocity_gradient
            IDENTITY_MAT22,                                 # deformation_gradient
            IDENTITY_MAT22,                                 # deformation_gradient_bar
            IDENTITY_MAT22,                                 # relative_deformation_gradient
            0.0,                                            # prev_jacobian_bar
            0.0,                                            # mass
            0.0,                                            # volume
            0.0,                                            # volume_init
            0.0,                                            # density
            0.0,                                            # pressure
            0.0,                                            # parameter_1
            0.0,                                            # parameter_2
            1,                                              # bulk_id
            zeros(UInt32, 4),                               # connected_grid_indices
        )
    end
end

