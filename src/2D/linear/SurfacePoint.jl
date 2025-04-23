
mutable struct SurfacePoint
    # Properties
    position::Vec2{Float64}
    surface_tension::Float64
    area::Float64
    area_init::Float64
    normal::Vec2{Float64}
    normal_init::Vec2{Float64}
    traction::Vec2{Float64}

    deformation_gradient::Mat22{Float64}

    # connectivity
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 8

    function SurfacePoint()
        new(ZERO_VEC2,                        # position
            0.0,                                            # surface_tension
            0.0,                                            # area
            0.0,                                            # area_init
            ZERO_VEC2,                        # normal
            ZERO_VEC2,                        # normal_init
            ZERO_VEC2,
            IDENTITY_MAT22,
            zeros(UInt32, 4),                               # connected_grid_indices
            )
    end
end
