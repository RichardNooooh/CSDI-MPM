
# Map the membrane points to the grid
function membraneToGrid!(grid::Grid, membrane_domain::Vector{MembranePoint}, time_scale_factor::Float64)::Nothing
    @threads for membrane_point::MembranePoint in membrane_domain
        for connected_grid_index::UInt8 in 1:membrane_point.connected_grid_array_length
            adj_grid_index::UInt32 = membrane_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, surface_gradi::Vec3{Float64} = surfCPDITri3ShapeFunctionAndGradient(membrane_point.vertices,
                adj_grid_point.position, grid.cell_length)

            I::Mat33{Float64} = IDENTITY_MAT33
            a::Float64, n::Vec3{Float64} = getTri3AreaAndNormal(membrane_point.vertices)

            grid_mass_term::Float64 = shape_value * membrane_point.mass
            grid_momentum_term::Vec3{Float64} = (shape_value * membrane_point.mass) .* membrane_point.velocity
            grid_force_term::Vec3{Float64} = time_scale_factor * membrane_point.surface_tension * ((I - n * n') * surface_gradi) * a
            lock(adj_grid_point.lock) do # * Note: it's possible to add gravity and traction to this! See paper.
                adj_grid_point.mass += grid_mass_term
                adj_grid_point.momentum += grid_momentum_term
                adj_grid_point.force -= grid_force_term
            end
        end
    end

    return nothing
end

# Update the momenta of the grid
function updateGrid!(grid_points::Vector{GridPoint}, dt::Float64)::Nothing
    @threads for grid_point::GridPoint in grid_points
        grid_point.momentum_next = grid_point.momentum + grid_point.force * dt

        # Fix stationary dirichlet boundary grid points
        if (grid_point.is_fixed[1] == true)
            grid_point.momentum = Vec3{Float64}(0.0, grid_point.momentum[2], grid_point.momentum[3])
            grid_point.momentum_next = Vec3{Float64}(0.0, grid_point.momentum_next[2], grid_point.momentum_next[3])
        end
        if (grid_point.is_fixed[2] == true)
            grid_point.momentum = Vec3{Float64}(grid_point.momentum[1], 0.0, grid_point.momentum[3])
            grid_point.momentum_next = Vec3{Float64}(grid_point.momentum_next[1], 0.0, grid_point.momentum_next[3])
        end
        if (grid_point.is_fixed[3] == true)
            grid_point.momentum = Vec3{Float64}(grid_point.momentum[1], grid_point.momentum[2], 0.0)
            grid_point.momentum_next = Vec3{Float64}(grid_point.momentum_next[1], grid_point.momentum_next[2], 0.0)
        end

        # Precalculate grid velocity for next step
        if grid_point.mass > EPSILON
            grid_point.velocity = grid_point.momentum ./ grid_point.mass
            grid_point.velocity_next = grid_point.momentum_next ./ grid_point.mass
        end
    end

    return nothing
end

function gridToMembrane!(membrane_domain::Vector{MembranePoint}, grid::Grid, alpha::Float64)::Nothing
    @threads for membrane_point::MembranePoint in membrane_domain
        grid_velocity_sum::Vec3{Float64} = ZERO_VEC3
        grid_velocity_increment_sum::Vec3{Float64} = ZERO_VEC3
        for connected_grid_index::UInt8 in 1:membrane_point.connected_grid_array_length
            adj_grid_index::UInt32 = membrane_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, surface_gradi::Vec3{Float64} = surfCPDITri3ShapeFunctionAndGradient(membrane_point.vertices,
                adj_grid_point.position, grid.cell_length)
            grid_velocity_sum += shape_value .* adj_grid_point.velocity_next
            grid_velocity_increment_sum += shape_value .* (adj_grid_point.velocity_next - adj_grid_point.velocity)
        end
        v::Vec3{Float64} = membrane_point.velocity
        membrane_point.velocity = alpha * (v + grid_velocity_increment_sum) + (1 - alpha) * grid_velocity_sum
    end

    return nothing
end

function updateMembraneVertices!(membrane_domain::Vector{MembranePoint}, grid::Grid, dt::Float64)::Nothing
    @threads for membrane_point::MembranePoint in membrane_domain
        for vertex_index::Int64 in 1:length(membrane_point.vertices)
            corner_position::Vec3{Float64} = membrane_point.vertices[vertex_index]
            corner_increment::Vec3{Float64} = ZERO_VEC3

            corner_adjacencies::Vec8{UInt32} = findConnectedGrid(corner_position, grid.cell_length, grid.num_nodes)
            for grid_index::UInt32 in corner_adjacencies
                grid_point::GridPoint = grid.points[grid_index]
                shape_value::Float64 = linearFunction(corner_position, grid_point.position, grid.cell_length)

                corner_increment += shape_value .* grid_point.velocity_next
            end
            membrane_point.vertices[vertex_index] += dt * corner_increment
        end
        membrane_point.connected_grid_array_length = getAllConnectedGrid!(membrane_point.connected_grid_indices,
            membrane_point.vertices,
            grid.cell_length, grid.num_nodes)

        membrane_point.area, membrane_point.normal = getTri3AreaAndNormal(membrane_point.vertices)
    end

    return nothing
end


function updateRigidGridFixed!(fixed_rigid_grid::Vector{GridPoint})::Nothing # this should be very small - should not use @threads here
    for fixed_grid_point::GridPoint in fixed_rigid_grid
        fixed_grid_point.momentum      = ZERO_VEC3
        # fixed_grid_point.momentum_next = ZERO_VEC3
        fixed_grid_point.force         = ZERO_VEC3
    end
end

function updateRigidVertices!(rigid_domain::Vector{RigidPoint}, grid::Grid, dt::Float64)::Nothing
    @threads for rigid_point::RigidPoint in rigid_domain
        velocity::Vec3{Float64} = rigid_point.velocity
        for vertex_index::Int64 in 1:length(rigid_point.vertices)
            rigid_point.vertices[vertex_index] += velocity * dt
        end

        rigid_point.connected_grid_array_length = getAllConnectedGrid!(rigid_point.connected_grid_indices, rigid_point.vertices, 
                                                                        grid.cell_length, grid.num_nodes)
    end

    return nothing
end


