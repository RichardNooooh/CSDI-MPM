
function surfaceToGrid!(grid::Grid, surface_domain::Vector{SurfacePoint}, t::Float64, t_load::Float64)::Nothing
    time_scale_factor::Float64 = t_load > 0.0 ? minimum((t_load, t)) / t_load : 1.0
    @threads for surface_point::SurfacePoint in surface_domain
        for connected_grid_index::UInt8 in 1:4
            adj_grid_index::UInt32 = surface_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            value::Float64, gradient::Vec2{Float64} = linearFunctionAndGradient(surface_point.position,
                adj_grid_point.position, grid.cell_length)

            I::Mat22{Float64} = IDENTITY_MAT22
            a::Float64, n::Vec2{Float64} = surface_point.area, surface_point.normal
            A::Float64 = surface_point.area_init

            surface_tension::Vec2{Float64} = time_scale_factor * surface_point.surface_tension * ((I - n * n') * gradient) * a
            traction::Vec2{Float64} = time_scale_factor * value * surface_point.traction * A
            grid_force_term::Vec2{Float64} = traction - surface_tension
            lock(adj_grid_point.lock) do
                adj_grid_point.force += grid_force_term
            end
        end
    end

    return nothing
end

# Map the grid to the material points
function gridToSurface!(surface_domain::Vector{SurfacePoint}, grid::Grid, dt::Float64)::Nothing
    @threads for surface_point::SurfacePoint in surface_domain
        velocity_gradient::Mat22{Float64} = ZERO_MAT22
        for connected_grid_index::UInt8 in 1:4
            adj_grid_index::UInt32 = surface_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_gradi::Vec2{Float64} = linearGradient(surface_point.position,
                adj_grid_point.position, grid.cell_length)

            grid_velocity::Vec2{Float64} = ZERO_VEC2
            if adj_grid_point.mass > EPSILON
                grid_velocity = adj_grid_point.momentum ./ adj_grid_point.mass
            end
            velocity_gradient += grid_velocity * shape_gradi'
        end
        I::Mat22{Float64} = IDENTITY_MAT22
        delF::Mat22{Float64} = I + (velocity_gradient .* dt)
        F::Mat22{Float64} = surface_point.deformation_gradient

        surface_point.deformation_gradient = delF * F
    end

    return nothing
end


# update corner positions
function updateSurfacePositions!(surface_domain::Vector{SurfacePoint}, grid::Grid, dt::Float64)::Nothing
    @threads for surface_point::SurfacePoint in surface_domain
        position_increment::Vec2{Float64} = ZERO_VEC2

        for grid_index::UInt32 in surface_point.connected_grid_indices
            grid_point::GridPoint = grid.points[grid_index]
            shape_value::Float64 = linearFunction(surface_point.position, grid_point.position, grid.cell_length)

            grid_velocity::Vec2{Float64} = ZERO_VEC2
            if grid_point.mass > EPSILON
                grid_velocity = grid_point.momentum ./ grid_point.mass
            end

            position_increment += (shape_value * dt) .* grid_velocity
        end
        surface_point.position += position_increment
        surface_point.connected_grid_indices = findConnectedGrid(surface_point.position,
            grid.cell_length, grid.num_nodes)

        # Nanson's formula
        F::Mat22{Float64} = surface_point.deformation_gradient
        J::Float64 = det(F)
        A::Float64, N::Vec2{Float64} = surface_point.area_init, surface_point.normal_init
        areanorm::Vec2{Float64} = J * A * inv(F)' * N

        a::Float64 = norm(areanorm)
        n::Vec2{Float64} = areanorm / a

        surface_point.area, surface_point.normal = a, n
    end

    return nothing
end
