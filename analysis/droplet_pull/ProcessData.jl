
function getLeftRightDroplets(material_domain::Vector{MaterialPoint},
    surface_domain::Vector{SurfacePoint},
    axis::Float64)::Tuple{Vector{MaterialPoint},Vector{MaterialPoint},Vector{SurfacePoint},Vector{SurfacePoint}}
    # creates vector of pointers to the material points on the left and right droplet domains

    left_material_domain::Vector{MaterialPoint} = Vector{MaterialPoint}(undef, 0)
    right_material_domain::Vector{MaterialPoint} = Vector{MaterialPoint}(undef, 0)

    for material_point::MaterialPoint in material_domain
        if abs(material_point.vertices[1][1]) < axis
            push!(left_material_domain, material_point)
        else
            push!(right_material_domain, material_point)
        end
    end

    @info("Found $(length(left_material_domain)) MPs in left, $(length(right_material_domain)) MPs in right droplets")

    left_surface_domain::Vector{SurfacePoint} = Vector{SurfacePoint}(undef, 0)
    right_surface_domain::Vector{SurfacePoint} = Vector{SurfacePoint}(undef, 0)

    for surface_point::SurfacePoint in surface_domain
        if abs(surface_point.vertices[1][1]) < axis
            push!(left_surface_domain, surface_point)
        else
            push!(right_surface_domain, surface_point)
        end
    end

    @info("Found $(length(left_surface_domain)) SPs in left, $(length(right_surface_domain)) SPs in right droplets")

    return left_material_domain, right_material_domain, left_surface_domain, right_surface_domain
end

function getTestGridpointsOnSymmetry(grid::Grid, axis::Float64)::Dict{UInt32,Tuple{GridPoint,GridPoint}}
    # creates new set of grid points on same positions as `grid` on x=`axis`
    # also records the corresponding grid.points index for recording on VTKs
    symAxisGridPoints::Dict{UInt32,Tuple{GridPoint,GridPoint}} = Dict{UInt32,Tuple{GridPoint,GridPoint}}()
    for (i::Int64, grid_point::GridPoint) in enumerate(grid.points)
        if abs(grid_point.position[1] - axis) < EPSILON
            newLeftGridPoint::GridPoint = GridPoint(grid_point.position)
            newRightGridPoint::GridPoint = GridPoint(grid_point.position)

            symAxisGridPoints[i] = (newLeftGridPoint, newRightGridPoint)
        end
    end

    @info("Found $(length(symAxisGridPoints)) grid points on the axis. Allocated $(Base.summarysize(symAxisGridPoints)) bytes.")

    return symAxisGridPoints
end

function oneSidedGridForceCalculation!(grid::Grid, isRight::Bool,
    material_domain_sided::Vector{MaterialPoint},
    surface_domain_sided::Vector{SurfacePoint},
    symAxisGridPoints::Dict{UInt32,Tuple{GridPoint,GridPoint}},
    time_scale_factor::Float64)::Nothing

    @threads for material_point::MaterialPoint in material_domain_sided
        for connected_grid_index::Int8 in 1:material_point.connected_grid_array_length
            adj_grid_index::UInt32 = material_point.connected_grid_indices[connected_grid_index]

            if !haskey(symAxisGridPoints, adj_grid_index)
                continue
            end

            axis_grid_point::GridPoint = symAxisGridPoints[adj_grid_index][1]
            if isRight
                axis_grid_point = symAxisGridPoints[adj_grid_index][2]
            end
            shape_value::Float64, shape_gradi::Vec2{Float64} = bulkCPDIShapeFunctionAndGradient(material_point.vertices,
                axis_grid_point.position, grid.cell_length)

            grid_mass_term::Float64 = shape_value * material_point.mass
            grid_momentum_term::Vec2{Float64} = (shape_value * material_point.mass) .* material_point.velocity
            grid_force_term::Vec2{Float64} = time_scale_factor * shape_value * material_point.force_external -
                                             material_point.volume .* (material_point.stress * shape_gradi)

            lock(axis_grid_point.lock) do
                axis_grid_point.mass += grid_mass_term
                axis_grid_point.momentum += grid_momentum_term
                axis_grid_point.force += grid_force_term
            end
        end
    end

    @threads for surface_point::SurfacePoint in surface_domain_sided
        for connected_grid_index::Int8 in 1:surface_point.connected_grid_array_length
            adj_grid_index::UInt32 = surface_point.connected_grid_indices[connected_grid_index]

            if !haskey(symAxisGridPoints, adj_grid_index)
                continue
            end

            axis_grid_point::GridPoint = symAxisGridPoints[adj_grid_index][1]
            if isRight
                axis_grid_point = symAxisGridPoints[adj_grid_index][2]
            end

            value::Float64, gradient::Vec2{Float64} = surfCPDIL2ShapeFunctionAndGradient(surface_point.vertices, 
                                                    axis_grid_point.position, grid.cell_length)
            
            I::Mat22{Float64} = IDENTITY_MAT22
            a::Float64, n::Vec2{Float64} = getL2AreaAndNormal(surface_point.vertices)
            A::Float64 = surface_point.area_init

            surface_tension::Vec2{Float64} = time_scale_factor * surface_point.surface_tension * ((I - n * n') * gradient) * a
            traction::Vec2{Float64} = time_scale_factor * value * surface_point.traction * A
            grid_force_term::Vec2{Float64} = traction - surface_tension
            lock(axis_grid_point.lock) do
                axis_grid_point.force += grid_force_term
            end
        end
    end
end


function writeVTKGridCustom!(grid::Grid,
    grid_data::GridVTKData, symAxisGridPoints::Dict{UInt32,Tuple{GridPoint,GridPoint}},
    plot_left_force::Array{Float64, 2}, plot_right_force::Array{Float64, 2},
    plot_force_over_time::Vector{Tuple{Float64, Vec2{Float64}, Vec2{Float64}}},
    pvd_grid::WriteVTK.CollectionFile,
    t_index::Int64, t_record::Int64, t::Float64,
    vtk_output_folder::String)::Nothing
    numGP::Int64 = length(grid.points)
    for gp_index = 1:numGP
        grid_point::GridPoint = grid.points[gp_index]

        # record all plotting information into these arrays
        grid_data.mass[gp_index] = grid_point.mass
        grid_data.momentum[gp_index, 1] = grid_point.momentum[1]
        grid_data.momentum[gp_index, 2] = grid_point.momentum[2]
        grid_data.force[gp_index, 1] = grid_point.force[1]
        grid_data.force[gp_index, 2] = grid_point.force[2]
    end
    fill!(plot_left_force, 0.0)
    fill!(plot_right_force, 0.0)
    total_left_force::Vec2{Float64} = ZERO_VEC2
    total_right_force::Vec2{Float64} = ZERO_VEC2
    for (gp_idx, (left_grid_point, right_grid_point)) in symAxisGridPoints
        plot_left_force[gp_idx, 1] = left_grid_point.force[1]
        plot_left_force[gp_idx, 2] = left_grid_point.force[2]

        plot_right_force[gp_idx, 1] = right_grid_point.force[1]
        plot_right_force[gp_idx, 2] = right_grid_point.force[2]

        total_left_force += left_grid_point.force
        total_right_force += right_grid_point.force
    end
    # println("left: $(total_left_force); right: $(total_right_force)")

    push!(plot_force_over_time, (t, total_left_force, total_right_force))

    # create vtk file
    vtk_grid(joinpath(vtk_output_folder, "grid_timestep_$(Int64(t_index / t_record))"),
        0.0:grid.cell_length[1]:grid.grid_length[1],
        0.0:grid.cell_length[2]:grid.grid_length[2]) do vtk
        vtk["Mass"] = grid_data.mass[:]
        vtk["Momentum"] = @views (grid_data.momentum[:, 1], grid_data.momentum[:, 2], grid_data.zero_cell)
        vtk["Force"] = @views (grid_data.force[:, 1], grid_data.force[:, 2], grid_data.zero_cell)
        vtk["Force From Left Droplet"] = @views (plot_left_force[:, 1], plot_left_force[:, 2], grid_data.zero_cell)
        vtk["Force From Right Droplet"] = @views (plot_right_force[:, 1], plot_right_force[:, 2], grid_data.zero_cell)
        pvd_grid[t] = vtk
    end

    return nothing
end

