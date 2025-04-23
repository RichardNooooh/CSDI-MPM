
# Map the material points to the grid
function materialToGrid!(grid::Grid, material_domain::Vector{MaterialPoint}, t::Float64, t_load::Float64)::Nothing
    time_scale_factor::Float64 = t_load > 0.0 ? minimum((t_load, t)) / t_load : 1.0
    @threads for material_point::MaterialPoint in material_domain
        for connected_grid_index::Int8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, grid.cell_length)

            grid_mass_term::Float64 = shape_value * material_point.mass
            grid_momentum_term::Vec2{Float64} = (shape_value * material_point.mass) .* material_point.velocity
            grid_force_term::Vec2{Float64} = time_scale_factor * shape_value * material_point.force_external -
                                             material_point.volume .* (material_point.stress * shape_gradi)

            lock(adj_grid_point.lock) do
                adj_grid_point.mass += grid_mass_term
                adj_grid_point.momentum += grid_momentum_term
                adj_grid_point.force += grid_force_term
            end
        end
    end

    return nothing
end

# Update the momenta of the grid
function updateGrid!(grid::Grid, dt::Float64)::Nothing
    @threads for grid_point::GridPoint in grid.points
        grid_point.momentum += grid_point.force * dt

        # Fix stationary dirichlet boundary grid points
        if (grid_point.is_fixed[1] == true)
            grid_point.momentum = Vec2{Float64}(0.0, grid_point.momentum[2])
            grid_point.force = Vec2{Float64}(0.0, grid_point.force[2])
        end
        if (grid_point.is_fixed[2] == true)
            grid_point.momentum = Vec2{Float64}(grid_point.momentum[1], 0.0)
            grid_point.force = Vec2{Float64}(grid_point.force[1], 0.0)
        end
    end

    return nothing
end

# Map the grid to the material points
function gridToMaterial!(material_domain::Vector{MaterialPoint}, grid::Grid, dt::Float64, applyConstitutiveEquation!::Function)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        velocity_gradient::Mat22{Float64} = ZERO_MAT22
        for connected_grid_index::UInt8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, grid.cell_length)

            # TODO this should probably be moved to gridUpdate()...
            grid_velocity::Vec2{Float64} = ZERO_VEC2
            grid_acceleration::Vec2{Float64} = ZERO_VEC2
            if adj_grid_point.mass > EPSILON
                grid_velocity = adj_grid_point.momentum ./ adj_grid_point.mass
                grid_acceleration = adj_grid_point.force ./ adj_grid_point.mass
            end

            material_point.velocity += (shape_value * dt) .* grid_acceleration

            velocity_gradient += grid_velocity * shape_gradi'
        end
        I::Mat22{Float64} = IDENTITY_MAT22
        material_point.velocity_gradient = velocity_gradient
        material_point.relative_deformation_gradient = I + (velocity_gradient .* dt)
        material_point.strain += velocity_gradient .* dt - I

        F::Mat22{Float64} = material_point.deformation_gradient
        material_point.deformation_gradient = material_point.relative_deformation_gradient * F

        applyConstitutiveEquation!(material_point)
    end

    return nothing
end

function RMohanFBarGridToMaterial!(material_domain::Vector{MaterialPoint}, grid::Grid, projection_grid::ProjectionGrid, dt::Float64, applyConstitutiveEquation!::Function)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        # Obtain the velocity gradient
        velocity_gradient::Mat22{Float64} = ZERO_MAT22
        for connected_grid_index::UInt8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, grid.cell_length)
            grid_velocity::Vec2{Float64} = ZERO_VEC2
            grid_acceleration::Vec2{Float64} = ZERO_VEC2
            if adj_grid_point.mass > EPSILON
                grid_velocity = adj_grid_point.momentum ./ adj_grid_point.mass
                grid_acceleration = adj_grid_point.force ./ adj_grid_point.mass
            end

            material_point.velocity += (shape_value * dt) .* grid_acceleration

            velocity_gradient += grid_velocity * shape_gradi'
        end
        material_point.velocity_gradient = velocity_gradient
        material_point.strain += velocity_gradient .* dt - I

        div_velocity::Float64 = tr(velocity_gradient)
        material_point.div_velocity = div_velocity

        # project to the grid
        # @assert grid.cell_length[1] == projection_grid.cell_length[1]

        material_point.connected_grid_indices = findConnectedGrid(material_point.position, material_point.side_lengths,
            projection_grid.cell_length, projection_grid.num_nodes)

        V::Float64 = material_point.volume
        for connected_grid_index::Int8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::ProjectionPoint = projection_grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, grid.cell_length)

            grid_volume_term::Float64 = shape_value * V
            avg_div_vel::Float64 = shape_value * div_velocity * V
            lock(adj_grid_point.lock) do
                adj_grid_point.avg_volume[material_point.bulk_ID] += grid_volume_term
                adj_grid_point.avg_div_vel[material_point.bulk_ID] += avg_div_vel
            end
        end
    end

    # Since we now obtained div_vel at the nodes, project it back
    @threads for material_point::MaterialPoint in material_domain
        projected_div_vel::Float64 = 0.0
        for connected_grid_index::Int8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::ProjectionPoint = projection_grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, grid.cell_length)

            avg_volume::Float64 = adj_grid_point.avg_volume[material_point.bulk_ID]
            avg_div_vel::Float64 = adj_grid_point.avg_div_vel[material_point.bulk_ID]
            if avg_volume > EPSILON
                projected_div_vel += shape_value * avg_div_vel / avg_volume
            end
        end
        ∇v::Mat22{Float64} = material_point.velocity_gradient
        div_v::Float64 = material_point.div_velocity
        I::Mat22{Float64} = IDENTITY_MAT22
        gradient_velocity_bar::Mat22{Float64} = ∇v + (1 / 3) * (projected_div_vel - div_v) * I

        material_point.velocity_gradient = gradient_velocity_bar

        Fbar::Mat22{Float64} = material_point.deformation_gradient_bar
        material_point.deformation_gradient_bar = (I + gradient_velocity_bar .* dt) * Fbar
        applyConstitutiveEquation!(material_point)
    end
end

function averageVolumesFBar!(material_domain::Vector{MaterialPoint}, projection_grid::ProjectionGrid, applyConstitutiveEquation!::Function)::Nothing
    # * MODIFIED F Bar Technique from Zhao et al (2022) - "Circumventing volumetric locking..."
    @threads for material_point::MaterialPoint in material_domain
        Jp::Float64 = det(material_point.deformation_gradient)
        V::Float64 = material_point.volume

        # TODO add connected_grid_indices for the projection grid
        # obtain connected grid points for projection grid
        material_point.connected_grid_indices = findConnectedGrid(material_point.position, material_point.side_lengths,
            projection_grid.cell_length, projection_grid.num_nodes)

        for connected_grid_index::Int8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::ProjectionPoint = projection_grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, projection_grid.cell_length)

            grid_volume_term::Float64 = shape_value * V
            grid_avg_jac::Float64 = shape_value * Jp * V
            lock(adj_grid_point.lock) do
                adj_grid_point.avg_volume[material_point.bulk_ID] += grid_volume_term
                adj_grid_point.avg_jacobian[material_point.bulk_ID] += grid_avg_jac
            end
        end
    end

    @threads for material_point::MaterialPoint in material_domain
        JBar::Float64 = 0.0
        for connected_grid_index::Int8 in 1:4
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::ProjectionPoint = projection_grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                adj_grid_point.position, material_point.side_lengths, projection_grid.cell_length)

            avg_volume::Float64 = adj_grid_point.avg_volume[material_point.bulk_ID]
            avg_jacobian::Float64 = adj_grid_point.avg_jacobian[material_point.bulk_ID]
            if avg_volume > EPSILON
                JBar += shape_value * avg_jacobian / avg_volume
            end
        end

        F::Mat22{Float64} = material_point.deformation_gradient
        J::Float64 = det(F)

        material_point.deformation_gradient_bar = (JBar / J)^(1 / 3) * F
        applyConstitutiveEquation!(material_point)
    end
end

# function avgCellFBar!(material_domain::Vector{MaterialPoint}, projection_grid::ProjectionGrid, applyConstitutiveEquation!::Function)::Nothing
#     @threads for material_point::MaterialPoint in material_domain
#         Jp::Float64 = det(material_point.deformation_gradient)

#         cell_index::UInt32 = getCellIndex(material_point.position, projection_grid.num_nodes, projection_grid.cell_length)
#         projection_point::ProjectionPoint = projection_grid.points[cell_index]

#         lock(projection_point.lock) do
#             projection_point.avg_volume[material_point.bulk_ID] += material_point.volume
#             projection_point.avg_jacobian[material_point.bulk_ID] += material_point.volume * Jp
#         end
#     end

#     @threads for material_point::MaterialPoint in material_domain
#         F::Mat22{Float64} = material_point.deformation_gradient
#         Jp = det(F)

#         cell_index::UInt32 = getCellIndex(material_point.position, projection_grid.num_nodes, projection_grid.cell_length)
#         projection_point::ProjectionPoint = projection_grid.points[cell_index]

#         J0::Float64 = projection_point.avg_jacobian[material_point.bulk_ID] / projection_point.avg_volume[material_point.bulk_ID]

#         material_point.deformation_gradient_bar = (J0 / Jp)^(1/3) * F
#         applyConstitutiveEquation!(material_point)
#     end
# end

# update corner positions
function updateMaterialPositions!(material_domain::Vector{MaterialPoint}, grid::Grid, dt::Float64, update_connected_grid::Bool)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        # if update_connected_grid
        #     material_point.connected_grid_indices = findConnectedGrid(material_point.position, material_point.side_lengths,
        #                                                                 grid.cell_length, grid.num_nodes)
        # end
        position_increment::Vec2{Float64} = ZERO_VEC2

        for grid_index::UInt32 in material_point.connected_grid_indices
            grid_point::GridPoint = grid.points[grid_index]
            shape_value::Float64, shape_gradi::Vec2{Float64} = GIMPShapeFunctionAndGradient(material_point.position,
                grid_point.position, material_point.side_lengths, grid.cell_length)
            grid_velocity::Vec2{Float64} = ZERO_VEC2
            if grid_point.mass > EPSILON
                grid_velocity = grid_point.momentum ./ grid_point.mass
            end

            position_increment += (shape_value * dt) .* grid_velocity
        end
        material_point.position += position_increment


        F::Mat22{Float64} = material_point.deformation_gradient
        J::Float64 = det(F)
        V_0::Float64 = material_point.volume_init
        material_point.volume = J * V_0
        # material_point.volume = getQ4Area(material_point.vertices)
        material_point.density = material_point.mass / material_point.volume

        # update cpgimp shape
        l_p::Vec2{Float64} = material_point.side_lengths_init
        material_point.side_lengths = Vec2{Float64}(l_p[1] * F[1, 1], l_p[2] * F[2, 2])

        material_point.connected_grid_indices = findConnectedGrid(material_point.position, material_point.side_lengths,
            grid.cell_length, grid.num_nodes)
    end

    return nothing
end
