
mutable struct SurfacePoint
    # Properties
    surface_tension::Float64
    area::Float64
    area_init::Float64
    normal::Vec3{Float64}
    normal_init::Vec3{Float64}
    traction::Vec3{Float64}

    # VTK file writing
    node_connectivity_vtk::Vec3{Int64}

    # connectivity
    connected_grid_array_length::UInt8     # denotes the final valid index in `connected_grid_indices`
    connected_grid_indices::Vector{UInt32} # one-time allocated Vector of length 24

    # CPDI-Tri3 vectex positions
    vertices::Vector{Vec3{Float64}}

    function SurfacePoint()
        new(0.0,                  # surface_tension
            0.0,                  # area
            0.0,                  # area_init
            ZERO_VEC3,            # normal
            ZERO_VEC3,            # normal_init
            ZERO_VEC3,            # traction
            Vec3{Int64}(0, 0, 0), # node_connectivity_vtk
            0,                    # connected_grid_array_length
            zeros(UInt32, 24),    # connected_grid_indices
            [ZERO_VEC3,           # vertices
             ZERO_VEC3,
             ZERO_VEC3])
    end
end

# utility for determining a surface point's position
# allows me to select surfaces without defining them in the mesh
function getSumOfVertexInDimension(surface_point::SurfacePoint, dim::Int64) # TODO refactor
    return sum((surface_point.vertices[1][dim], surface_point.vertices[2][dim], surface_point.vertices[3][dim]))
end
