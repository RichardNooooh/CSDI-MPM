
function trackFarthestCorner(material_domain::Vector{MaterialPoint})::Tuple{MaterialPoint, Int64}
    for material_point::MaterialPoint in material_domain
        for (index::Int64, vertex::Vec2{Float64}) in enumerate(material_point.vertices)
            if sum(abs.(vertex - Vec2{Float64}(1.0, 1.0))) < EPSILON
                return material_point, index
            end
        end
    end
end


function getElastocapillaryNumber(tracked_mp::MaterialPoint, tracked_vertex::Int64, arguments::Dict{String, Any})::Float64
    final_position::Vec2{Float64} = tracked_mp.vertices[tracked_vertex]
    r::Float64 = norm(final_position)

    l::Float64 = 2*r
    γ::Float64 = arguments["surface_tension"] * arguments["surface_tension_factor"]
    E::Float64 = arguments["elastic_modulus"]

    return γ / (E * l)
end

