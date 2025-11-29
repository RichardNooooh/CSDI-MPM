
"""
Calculates first and second Lamé parameters given the elastic modulus and Poisson ratio.

# Arguments
- `E::Float64`: elastic modulus
- `v::Float64`: Poisson ratio

# Output
- Returns tuple of first and second Lamé Parameters
"""
function getLaméParameters(E::Float64, v::Float64)::Tuple{Float64, Float64}
    λ::Float64 = E * v / ((1 + v) * (1 - 2*v))  # 1st Lamé Parameter
    μ::Float64 = 0.5 * E / (1 + v)              # shear modulus / 2nd Lamé Parameter

    return λ, μ
end


"""
Modifies stress tensor with an hyperelastic constitutive equation. Assumes
`parameter_1` and `parameter_2` are the 1st and 2nd Lame parameters.

# Arguments
- `material_point::MaterialPoint`: MaterialPoint object to modify

# Output
- Nothing
"""
function hyperelasticConstitutiveEquation!(material_point::MaterialPoint)::Nothing
    I::Mat33{Float64} = IDENTITY_MAT33
    F::Mat33{Float64} = material_point.deformation_gradient
    J::Float64 = det(F)
    λ::Float64 = material_point.parameter_1
    μ::Float64 = material_point.parameter_2

    p::Float64 = -λ / 2 * (J - 1/J)
    material_point.stress = -p * I +  μ * J^(-5/3) * (F * F' - (1/3) * (tr(F' * F)) * I)
    material_point.pressure = p
    return nothing
end


"""
Modifies stress tensor with an viscous fluid constitutive equation. Assumes
`parameter_1` is the bulk modulus and `parameter_2` is the dynamic viscosity.

# Arguments
- `material_point::MaterialPoint`: MaterialPoint object to modify

# Output
- Nothing
"""
function fluidConstitutiveEquation!(material_point::MaterialPoint)::Nothing
    I::Mat33{Float64} = IDENTITY_MAT33
    F::Mat33{Float64} = material_point.deformation_gradient
    J::Float64 = det(F)

    K::Float64 = material_point.parameter_1
    v::Float64 = material_point.parameter_2
    p::Float64 = -K / 2 * (J - 1 / J)     # hydrostatic pressure

    L::Mat33{Float64} = material_point.velocity_gradient
    e::Mat33{Float64} = 0.5 * (L + L')                              # strain rate tensor: 1/2 L + L^T
    d::Mat33{Float64} = e - (1 / 3) * tr(e) * I                     # deviatoric part of the strain rate

    material_point.stress = -p * I + 2 * v * d
    material_point.pressure = p

    return nothing
end

