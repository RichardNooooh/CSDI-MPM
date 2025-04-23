
# Map the material points to the grid
function materialToGrid!(grid::Grid, material_domain::Vector{MaterialPoint}, time_scale_factor::Float64)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        for connected_grid_index::Int8 in 1:material_point.connected_grid_array_length
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec3{Float64} = bulkCPDIShapeFunctionAndGradient(material_point.vertices, 
                                                                    adj_grid_point.position, grid.cell_length)
            
            grid_mass_term::Float64 = shape_value * material_point.mass
            grid_momentum_term::Vec3{Float64} = (shape_value * material_point.mass) .* material_point.velocity
            grid_force_term::Vec3{Float64} = time_scale_factor * shape_value * material_point.force_external - 
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

# Map the surface points to the grid
function surfaceToGrid!(grid::Grid, surface_domain::Vector{SurfacePoint}, time_scale_factor::Float64)::Nothing
    @threads for surface_point::SurfacePoint in surface_domain
        for connected_grid_index::UInt8 in 1:surface_point.connected_grid_array_length
            adj_grid_index::UInt32 = surface_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            value::Float64, gradient::Vec3{Float64} = surfCPDITri3ShapeFunctionAndGradient(surface_point.vertices, 
                                                    adj_grid_point.position, grid.cell_length)
            
            I::Mat33{Float64} = IDENTITY_MAT33
            a::Float64, n::Vec3{Float64} = getTri3AreaAndNormal(surface_point.vertices)
            A::Float64 = surface_point.area_init

            surface_tension::Vec3{Float64} = time_scale_factor * surface_point.surface_tension * ((I - n * n') * gradient) * a
            traction::Vec3{Float64} = time_scale_factor * value * surface_point.traction * A
            grid_force_term::Vec3{Float64} = traction - surface_tension
            lock(adj_grid_point.lock) do
                adj_grid_point.force += grid_force_term
            end
        end
    end
    
    return nothing
end

# Update the momenta of the grid
function updateGrid!(grid::Grid, dt::Float64)::Nothing
    @threads for grid_point::GridPoint in grid.points
        grid_point.momentum_next = grid_point.momentum + grid_point.force * dt

        # Fix stationary dirichlet boundary grid points
        if (grid_point.is_fixed[1] == true)
            grid_point.momentum      = Vec3{Float64}(0.0                        , grid_point.momentum[2]        , grid_point.momentum[3])
            grid_point.momentum_next = Vec3{Float64}(0.0                        , grid_point.momentum_next[2]   , grid_point.momentum_next[3])
        end
        if (grid_point.is_fixed[2] == true)
            grid_point.momentum      = Vec3{Float64}(grid_point.momentum[1]     , 0.0                           , grid_point.momentum[3])
            grid_point.momentum_next = Vec3{Float64}(grid_point.momentum_next[1], 0.0                           , grid_point.momentum_next[3])
        end
        if (grid_point.is_fixed[3] == true)
            grid_point.momentum      = Vec3{Float64}(grid_point.momentum[1]     , grid_point.momentum[2]        , 0.0)
            grid_point.momentum_next = Vec3{Float64}(grid_point.momentum_next[1], grid_point.momentum_next[2]   , 0.0)
        end

        # Precalculate grid velocity for next step
        if grid_point.mass > EPSILON
            grid_point.velocity = grid_point.momentum ./ grid_point.mass
            grid_point.velocity_next = grid_point.momentum_next ./ grid_point.mass
        end
    end

    return nothing
end

# Map the grid to the material points
function gridToMaterial!(material_domain::Vector{MaterialPoint}, grid::Grid, dt::Float64, alpha::Float64, applyConstitutiveEquation!::Function)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        grid_velocity_sum::Vec3{Float64} = ZERO_VEC3
        grid_velocity_increment_sum::Vec3{Float64} = ZERO_VEC3
        velocity_gradient::Mat33{Float64} = ZERO_MAT33
        for connected_grid_index::UInt8 in 1:material_point.connected_grid_array_length
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            shape_value::Float64, shape_gradi::Vec3{Float64} = bulkCPDIShapeFunctionAndGradient(material_point.vertices, 
                                                                    adj_grid_point.position, grid.cell_length)

            grid_velocity_sum += shape_value .* adj_grid_point.velocity_next
            grid_velocity_increment_sum += shape_value .* (adj_grid_point.velocity_next - adj_grid_point.velocity)

            velocity_gradient += adj_grid_point.velocity_next * shape_gradi'
        end
        v::Vec3{Float64} = material_point.velocity
        material_point.velocity = alpha * (v + grid_velocity_increment_sum) + (1-alpha) * grid_velocity_sum
        
        I::Mat33{Float64} = IDENTITY_MAT33
        material_point.velocity_gradient = velocity_gradient
        material_point.relative_deformation_gradient = I + (velocity_gradient .* dt)
        
        F::Mat33{Float64} = material_point.deformation_gradient
        material_point.deformation_gradient = material_point.relative_deformation_gradient * F
        applyConstitutiveEquation!(material_point)
    end

    return nothing
end

# update corner positions
function updateMaterialVertices!(material_domain::Vector{MaterialPoint}, grid::Grid, dt::Float64)::Nothing
    @threads for material_point::MaterialPoint in material_domain
        for vertex_index::Int64 in 1:length(material_point.vertices)
            corner_position::Vec3{Float64} = material_point.vertices[vertex_index]
            corner_increment::Vec3{Float64} = ZERO_VEC3

            corner_adjacencies::Vec8{UInt32} = findConnectedGrid(corner_position, grid.cell_length, grid.num_nodes)
            for grid_index::UInt32 in corner_adjacencies
                grid_point::GridPoint = grid.points[grid_index]
                shape_value::Float64 = linearFunction(corner_position, grid_point.position, grid.cell_length)

                corner_increment += (shape_value * dt) .* grid_point.velocity_next
            end
            material_point.vertices[vertex_index] += corner_increment
        end
        material_point.connected_grid_array_length = getAllConnectedGrid!(material_point.connected_grid_indices, 
                                                                                material_point.vertices, 
                                                                                grid.cell_length, grid.num_nodes)
        
        F::Mat33{Float64} = material_point.deformation_gradient
        J::Float64 = det(F)
        V_0::Float64 = material_point.volume_init
        material_point.volume = J * V_0
        material_point.density = material_point.mass / material_point.volume
    end

    return nothing
end

function updateSurfaceVertices!(surface_domain::Vector{SurfacePoint}, grid::Grid, dt::Float64)::Nothing
    @threads for surface_point::SurfacePoint in surface_domain
        for vertex_index::Int64 in 1:length(surface_point.vertices)
            corner_position::Vec3{Float64} = surface_point.vertices[vertex_index]
            corner_increment::Vec3{Float64} = ZERO_VEC3

            corner_adjacencies::Vec8{UInt32} = findConnectedGrid(corner_position, grid.cell_length, grid.num_nodes)
            for grid_index::UInt32 in corner_adjacencies
                grid_point::GridPoint = grid.points[grid_index]
                shape_value::Float64 = linearFunction(corner_position, grid_point.position, grid.cell_length)

                corner_increment += (shape_value * dt) .* grid_point.velocity_next
            end
            surface_point.vertices[vertex_index] += corner_increment
        end
        surface_point.connected_grid_array_length = getAllConnectedGrid!(surface_point.connected_grid_indices, 
                                                                                surface_point.vertices, 
                                                                                grid.cell_length, grid.num_nodes)
        
        surface_point.area, surface_point.normal = getTri3AreaAndNormal(surface_point.vertices)
    end

    return nothing
end


function updateRigidGrid!(grid::Grid, rigid_domain::Vector{RigidPoint})::Nothing
    @threads for rigid_point::RigidPoint in rigid_domain
        for connected_grid_index::Int8 in 1:rigid_point.connected_grid_array_length
            adj_grid_index::Int64 = rigid_point.connected_grid_indices[connected_grid_index]
            adj_grid_point::GridPoint = grid.points[adj_grid_index]

            adj_grid_point.momentum = adj_grid_point.mass * rigid_point.velocity
            adj_grid_point.force = ZERO_VEC2
        end
    end

    return nothing
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
