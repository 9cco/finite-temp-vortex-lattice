############################################################################################################################
#                               XY -> Chiral transformation // SubCuboidModule
#__________________________________________________________________________________________________________________________#
############################################################################################################################
# Functions for converting the xy-basis values of {u₁, u₂, θ₁, θ₂} to the Chiral basis values
# {u⁺, u⁻, θ⁺, θ⁻}.

# -----------------------------------------------------------------------------------------------------------
# Takes a number x and returns +1 if x >= 0, -1 if x < 0
function sgn(x::R) where R <: Real
    if x == 0
        return 1
    else
        return sign(x)
    end
end

# -----------------------------------------------------------------------------------------------------------
# Calculates the chiral amplitudes
function findu⁺(θ₁::R1, θ₂::R2, u₁::R3, u₂::R3) where {R1 <: Real, R2 <: Real, R3 <:Real}
    √((u₁^2 + u₂^2)/2 + u₁*u₂*sin(θ₁-θ₂))
end
function findu⁻(θ₁::R1, θ₂::R2, u₁::R3, u₂::R3) where {R1 <: Real, R2 <: Real, R3 <:Real}
    √((u₁^2 + u₂^2)/2 - u₁*u₂*sin(θ₁-θ₂))
end
function findu⁺(ϕ::LatticeSite)
    findu⁺(ϕ.θ⁺, ϕ.θ⁻, ϕ.u⁺, ϕ.u⁻)
end
function findu⁻(ϕ::LatticeSite)
    findu⁻(ϕ.θ⁺, ϕ.θ⁻, ϕ.u⁺, ϕ.u⁻)
end


# -----------------------------------------------------------------------------------------------------------
# Calculates θ⁺ when we assume that u₁ = u₂. When the phase difference θ₁-θ₂ is close to 3π/2 modulo 2π
# the chiral amplitude goes to zero along with the normal values used to calculate cos and sine creating
# numerical error. In this case a series expansion in ϵ = mod(θ₁-θ₂, 2π) - 3π/2 is used. The expresions
# for sine and cos are expanded to 4th order in ϵ.
function findθ⁺(θ₁::R1, θ₂::R2, u₁::R3, u₂::R3; ϵ = 1e-2) where {R1 <: Real, R2 <: Real, R3 <: Real}
    δ = mod(θ₁-θ₂,2π) - 3π/2
    
    if abs(δ) < ϵ
        arccosθ⁺ = acos(-sgn(δ)*sin(θ₁)*(1-δ^2/8+δ^4/384) + abs(δ)/2*cos(θ₁)*(1-δ^2/24))
        sinθ⁺ = sgn(δ)*cos(θ₁)*(1-δ^2/8) + abs(δ)/2*sin(θ₁)*(1-δ^2/24)
        if sinθ⁺ <= 0
            return -arccosθ⁺
        else
            return arccosθ⁺
        end
    else
        u⁺ = findu⁺(θ₁, θ₂, u₁, u₂)
        cosθ⁺_num = u₁*cos(θ₁) - u₂*sin(θ₂)
        sinθ⁺_num = u₁*sin(θ₁) + u₂*cos(θ₂)
        if cosθ⁺_num == 0
            return sgn(sinθ⁺_num)*π/2
        else
            tanθ⁺ = sinθ⁺_num/cosθ⁺_num
            arctanθ⁺ = atan(tanθ⁺)
            if cosθ⁺_num > 0
                return arctanθ⁺
            else
                return arctanθ⁺ - sgn(tanθ⁺)*π
            end
        end
    end
end
function findθ⁺(ϕ::LatticeSite)
    findθ⁺(ϕ.θ⁺, ϕ.θ⁻, ϕ.u⁺, ϕ.u⁻)
end

# Same as above, but in this case, the corresponding chiral amplitude u⁻ becomes zero when θ₁-θ₂ approaches
# π/2, thus the definition of the smallness parameter becomes ϵ = mod(θ₁-θ₂,2π) - π/2.
function findθ⁻(θ₁::R1, θ₂::R2, u₁::R3, u₂::R3; ϵ = 1e-2) where {R1 <: Real, R2 <: Real, R3 <: Real}
    δ = mod(θ₁-θ₂,2π) - π/2
    
    if abs(δ) < ϵ
        arccosθ⁻ = acos(-sgn(δ)*sin(θ₁)*(1-δ^2/8+δ^4/384) + abs(δ)/2*cos(θ₁)*(1-δ^2/24))
        sinθ⁻ = sgn(δ)*cos(θ₁)*(1-δ^2/8) + abs(δ)/2*sin(θ₁)*(1-δ^2/24)
        if sinθ⁻ <= 0
            return -arccosθ⁻
        else
            return arccosθ⁻
        end
    else
        u⁻ = findu⁻(θ₁, θ₂, u₁, u₂)
        cosθ⁻_num = u₁*cos(θ₁) + u₂*sin(θ₂)
        sinθ⁻_num = u₁*sin(θ₁) - u₂*cos(θ₂)
        if cosθ⁻_num == 0
            return sgn(sinθ⁻_num)*π/2
        else
            tanθ⁻ = sinθ⁻_num/cosθ⁻_num
            arctanθ⁻ = atan(tanθ⁻)
            if cosθ⁻_num > 0
                return arctanθ⁻
            else
                return arctanθ⁻ - sgn(tanθ⁻)*π
            end
        end
    end
end
function findθ⁻(ϕ::LatticeSite)
    findθ⁻(ϕ.θ⁺, ϕ.θ⁻, ϕ.u⁺, ϕ.u⁻)
end
