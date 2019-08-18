############################################################################################################################
#                               Functions for calculating Helicity Modulus
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# General helicity modulus for a twist θ± → θ± + δ(r⋅ω)a± in the ω direction
function twoCompHelMod(a⁺::R1, a⁻::R1, Υ₁₀::R2, Υ₀₁::R2, Υ₁₁::R2) where {R1<:Real, R2<:Real}
    (a⁺ - a⁻)(a⁺*Υ₁₀ - a⁻*Υ₀₁) + a⁺*a⁻*Υ₁₁
end

# -----------------------------------------------------------------------------------------------------------
# Helicity modulus definition, which is dependent on already having calculated thermal averages of the
# Hamiltonian derivatives
function helMod(β::R1, L₁::I, L₂::I, L₃::I, d²H_avg::R2, dH²_avg::R2, dH_avg::R2) where {R1<:Real, I<:Int, R2<:Real}
    (d²H_avg - β*(dH²_avg - dH_avg^2))/(L₁*L₂*L₃)
end

# -----------------------------------------------------------------------------------------------------------
# Calculates density of the first derivative of the free energy w.r.t. a twist in x direction.
# Since there are two phases the function calculates for a general twist coverned by a⁺ and a⁻
# First we have a function that converts from the xy-basis to the chiral basis values.
function firstDerivativeDensityX(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    firstDerivativeDensityX(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                           findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), A₁, A₂, x, y, ν, κ₅)
end
function firstDerivativeDensityX(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, A₁::R, A₂::R, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    fk_x = 2*(u⁺ᵣ₊₁*u⁺*sin(θ⁺ᵣ₊₁-θ⁺-A₁)a⁺ 
          + u⁻ᵣ₊₁*u⁻*sin(θ⁻ᵣ₊₁-θ⁻-A₁)a⁻)
    fmc_x = (ν+1)*(u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁ - θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))
               + u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁ - θ⁺ - A₁)*(a⁻-x*(a⁺-a⁻))
               - x*(a⁺-a⁻)*u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂ - θ⁻ - A₂)
               + x*(a⁺-a⁻)*u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂ - θ⁺ - A₂))
    fmg_x = (1-ν)*(u⁺ᵣ₊₁*u⁻ᵣ₊₂*cos(θ⁺ᵣ₊₁ - θ⁻ᵣ₊₂ - (A₁ - A₂))*(a⁺+x*(a⁺-a⁻)) 
                 + u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁ - θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻)) 
                 + x*(a⁺-a⁻)*(u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂ - θ⁻ - A₂) + u⁺*u⁻*cos(θ⁺-θ⁻))
                 + u⁻ᵣ₊₁*u⁺ᵣ₊₂*cos(θ⁻ᵣ₊₁ - θ⁺ᵣ₊₂ - (A₁ - A₂))*(-a⁻+x*(a⁺-a⁻)) 
                 + u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁ - θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻)) 
                 + x*(a⁺-a⁻)*(u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂ - θ⁺ - A₂) + u⁻*u⁺*cos(θ⁻-θ⁺)))
    fv_x = -ν*(u⁺*u⁻)^2*x*(a⁺-a⁻)*sin(2(θ⁺ - θ⁻))

    fk_x+fmc_x+fmg_x+fv_x
end
# Same but Y-direction
function firstDerivativeDensityY(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    firstDerivativeDensityY(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                           findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), A₁, A₂, x, y, ν, κ₅)
end

function firstDerivativeDensityY(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, A₁::R, A₂::R, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    fk_y = 2*(u⁺ᵣ₊₂*u⁺*sin(θ⁺ᵣ₊₂-θ⁺-A₂)*a⁺ 
              + u⁻ᵣ₊₂*u⁻*sin(θ⁻ᵣ₊₂-θ⁻-A₂)a⁻)
    fmc_y = -(ν+1)*(u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂ - θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))
                    + u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂ - θ⁺ - A₂)*(a⁻-y*(a⁺-a⁻))
               - y*(a⁺-a⁻)*u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁ - θ⁻ - A₁)
               + y*(a⁺-a⁻)*u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁ - θ⁺ - A₁))
    fmg_y = (1-ν)*(u⁺ᵣ₊₂*u⁻ᵣ₊₁*cos(θ⁺ᵣ₊₂ - θ⁻ᵣ₊₁ - (A₂ - A₁))*(-a⁺+y*(a⁺-a⁻)) 
                   + u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂ - θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻)) 
                   + y*(a⁺-a⁻)*(u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁ - θ⁻ - A₁) + u⁺*u⁻*cos(θ⁺-θ⁻))
                   + u⁻ᵣ₊₂*u⁺ᵣ₊₁*cos(θ⁻ᵣ₊₂ - θ⁺ᵣ₊₁ - (A₂ - A₁))*(a⁻+y*(a⁺-a⁻)) 
                   + u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂ - θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻)) 
                   + y*(a⁺-a⁻)*(u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁ - θ⁺ - A₁) + u⁻*u⁺*cos(θ⁻-θ⁺)))
    fv_y = -ν*(u⁺*u⁻)^2*y*(a⁺-a⁻)*sin(2(θ⁺ - θ⁻))

    fk_y+fmc_y+fmg_y+fv_y
end
# Same but Z-direction
function  firstDerivativeDensityZ(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, A::P, z::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₃ = ϕ.A₃
    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    firstDerivativeDensityZ(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                            findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), findu⁺(ϕᵣ₊₃), findu⁻(ϕᵣ₊₃), findθ⁺(ϕᵣ₊₃), findθ⁻(ϕᵣ₊₃),
                            A₁, A₂, A₃, z, ν, κ₅)
end

function firstDerivativeDensityZ(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, u⁺ᵣ₊₃::R, u⁻ᵣ₊₃::R, θ⁺ᵣ₊₃::R, θ⁻ᵣ₊₃::R, A₁::R, A₂::R, A₃::R, z::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    dk_z = 2*(u⁺ᵣ₊₃*u⁺*sin(θ⁺ᵣ₊₃-θ⁺-A₃)*a⁺ 
              + u⁻ᵣ₊₃*u⁻*sin(θ⁻ᵣ₊₃-θ⁻-A₃)a⁻)
    dmc_z = (ν+1)*z*((a⁺-a⁻)*(u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁ - θ⁻ - A₁) 
                              - u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂ - θ⁻ - A₂)) 
                   + (a⁻-a⁺)*(u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁ - θ⁺ - A₁) 
                              - u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂ - θ⁺ - A₂)))
    dmg_z = -(ν-1)z*(a⁺-a⁻)*(u⁺ᵣ₊₁*u⁻ᵣ₊₂*cos(θ⁺ᵣ₊₁-θ⁻ᵣ₊₂-(A₁-A₂)) + u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁-θ⁻-A₁) + u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂-θ⁻-A₂)
                             + u⁺*u⁻*cos(θ⁺-θ⁻)
                            u⁻ᵣ₊₁*u⁺ᵣ₊₂*cos(θ⁻ᵣ₊₁-θ⁺ᵣ₊₂-(A₁-A₂)) + u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁-θ⁺-A₁) + u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂-θ⁺-A₂)
                            + u⁻*u⁺*cos(θ⁻-θ⁺))
    dv_z = -ν*(u⁺*u⁻)^2*z*(a⁺-a⁻)*sin(2(θ⁺-θ⁻))

    dk_z+dmc_z+dmg_z+dv_z
end



# -----------------------------------------------------------------------------------------------------------
# Calculates the density of second derivative of the free energy w.r.t. a twist in x direction
# for a general (a⁺, a⁻) twist.
function secondDerivativeDensityX(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    secondDerivativeDensityX(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                            findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), A₁, A₂, x, ν, κ₅)
end

function secondDerivativeDensityX(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, A₁::R, A₂::R, x::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    ddk_x = 2*(u⁺ᵣ₊₁*u⁺*cos(θ⁺ᵣ₊₁-θ⁺-A₁)(a⁺)^2
               + u⁻ᵣ₊₁*u⁻*cos(θ⁻ᵣ₊₁-θ⁻-A₁)(a⁻)^2)
    ddmc_x = (ν+1)*(u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁ - θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))^2
               + u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁ - θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻))^2
               - x^2*(a⁺-a⁻)^2*u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂ - θ⁻ - A₂)
               - x^2*(a⁺-a⁻)^2*u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂ - θ⁺ - A₂))
    ddmg_x = (ν-1)*(u⁺ᵣ₊₁*u⁻ᵣ₊₂*sin(θ⁺ᵣ₊₁ - θ⁻ᵣ₊₂ - (A₁ - A₂))*(a⁺+x*(a⁺-a⁻))^2 
                 + u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁ - θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))^2 
                 + x^2*(a⁺-a⁻)^2*(u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂ - θ⁻ - A₂) + u⁺*u⁻*sin(θ⁺-θ⁻))
                 - u⁻ᵣ₊₁*u⁺ᵣ₊₂*sin(θ⁻ᵣ₊₁ - θ⁺ᵣ₊₂ - (A₁ - A₂))*(-a⁻+x*(a⁺-a⁻))^2 
                 - u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁ - θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻))^2 
                 - x^2*(a⁺-a⁻)^2*(u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂ - θ⁺ - A₂) + u⁻*u⁺*sin(θ⁻-θ⁺)))
    ddv_x = -ν*(u⁺*u⁻)^2*x^2*(a⁺-a⁻)^2*cos(2(θ⁺ - θ⁻))

    ddk_x+ddmc_x+ddmg_x+ddv_x
end
# Same but in Y-direction
function secondDerivativeDensityY(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    secondDerivativeDensityY(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                            findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), A₁, A₂, y, ν, κ₅)
end

function secondDerivativeDensityY(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, A₁::R, A₂::R, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    ddk_y = 2*(u⁺ᵣ₊₂*u⁺*cos(θ⁺ᵣ₊₂-θ⁺-A₂)(a⁺)^2
               + u⁻ᵣ₊₂*u⁻*cos(θ⁻ᵣ₊₂-θ⁻-A₂)(a⁻)^2)
    ddmc_y = -(ν+1)*(u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂ - θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))^2
                     + u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂ - θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻))^2
               - y^2*(a⁺-a⁻)^2*u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁ - θ⁻ - A₁)
               - y^2*(a⁺-a⁻)^2*u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁ - θ⁺ - A₁))
    ddmg_y = (ν-1)*(u⁺ᵣ₊₂*u⁻ᵣ₊₁*sin(θ⁺ᵣ₊₂ - θ⁻ᵣ₊₁ - (A₂ - A₁))*(a⁺-y*(a⁺-a⁻))^2 
                    + u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂ - θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))^2 
                 + y^2*(a⁺-a⁻)^2*(u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁ - θ⁻ - A₁) + u⁺*u⁻*sin(θ⁺-θ⁻))
                 - u⁻ᵣ₊₂*u⁺ᵣ₊₁*sin(θ⁻ᵣ₊₂ - θ⁺ᵣ₊₁ - (A₂ - A₁))*(-a⁻-y*(a⁺-a⁻))^2 
                 - u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂ - θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻))^2 
                 - y^2*(a⁺-a⁻)^2*(u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁ - θ⁺ - A₁) + u⁻*u⁺*sin(θ⁻-θ⁺)))
    ddv_y = -ν*(u⁺*u⁻)^2*y^2*(a⁺-a⁻)^2*cos(2(θ⁺ - θ⁻))

    ddk_y+ddmc_y+ddmg_y+ddv_y
end
# Same but in Z-direction
function secondDerivativeDensityZ(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, A::P, z::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁
    A₃ = ϕ.A₃

    secondDerivativeDensityZ(a⁺, a⁻, findu⁺(ϕ), findu⁻(ϕ), findθ⁺(ϕ), findθ⁻(ϕ), findu⁺(ϕᵣ₊₁), findu⁻(ϕᵣ₊₁), findθ⁺(ϕᵣ₊₁), findθ⁻(ϕᵣ₊₁),
                             findu⁺(ϕᵣ₊₂), findu⁻(ϕᵣ₊₂), findθ⁺(ϕᵣ₊₂), findθ⁻(ϕᵣ₊₂), findu⁺(ϕᵣ₊₃), findu⁻(ϕᵣ₊₃), findθ⁺(ϕᵣ₊₃), findθ⁻(ϕᵣ₊₃),
                             A₁, A₂, z, ν, κ₅)
end
function secondDerivativeDensityZ(a⁺::R, a⁻::R, u⁺::R, u⁻::R, θ⁺::R, θ⁻::R, u⁺ᵣ₊₁::R, u⁻ᵣ₊₁::R, θ⁺ᵣ₊₁::R, θ⁻ᵣ₊₁, u⁺ᵣ₊₂::R, u⁻ᵣ₊₂::R, θ⁺ᵣ₊₂::R
                                 θ⁻ᵣ₊₂::R, u⁺ᵣ₊₃::R, u⁻ᵣ₊₃::R, θ⁺ᵣ₊₃::R, θ⁻ᵣ₊₃::R, A₁::R, A₂::R, A₃::R, z::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    ddk_z = 2*(u⁺ᵣ₊₃*u⁺*cos(θ⁺ᵣ₊₃-θ⁺-A₃)*(a⁺)^2 + u⁻ᵣ₊₃*u⁻*cos(θ⁻ᵣ₊₃-θ⁻-A₃)*(a⁻)^2)
    ddmc_z = (ν+1)*z^2*(a⁺-a⁻)^2*(u⁺ᵣ₊₁*u⁻*cos(θ⁺ᵣ₊₁-θ⁻-A₁) - u⁺ᵣ₊₂*u⁻*cos(θ⁺ᵣ₊₂-θ⁻-A₂)
                                  + u⁻ᵣ₊₁*u⁺*cos(θ⁻ᵣ₊₁-θ⁺-A₁) - u⁻ᵣ₊₂*u⁺*cos(θ⁻ᵣ₊₂-θ⁺-A₂))
    ddmg_z = (ν-1)*z^2*(a⁺-a⁻)^2*(u⁺ᵣ₊₁*u⁻ᵣ₊₂*sin(θ⁺ᵣ₊₁-θ⁻ᵣ₊₂ - (A₁-A₂)) + u⁺ᵣ₊₁*u⁻*sin(θ⁺ᵣ₊₁-θ⁻-A₁) 
                                  + u⁺ᵣ₊₂*u⁻*sin(θ⁺ᵣ₊₂-θ⁻-A₂) + u⁺*u⁻*sin(θ⁺-θ⁻)
                                -(u⁻ᵣ₊₁*u⁺ᵣ₊₂*sin(θ⁻ᵣ₊₁-θ⁺ᵣ₊₂ - (A₁-A₂)) + u⁻ᵣ₊₁*u⁺*sin(θ⁻ᵣ₊₁-θ⁺-A₁) 
                                  + u⁻ᵣ₊₂*u⁺*sin(θ⁻ᵣ₊₂-θ⁺-A₂) + u⁻*u⁺*sin(θ⁻-θ⁺)))
    ddv_z = -ν*z^2*(a⁺-a⁻)^2*(u⁺*u⁻)^2*cos(2(θ⁺-θ⁻))

    ddk_z+ddmc_z+ddmg_z+ddv_z
end


# -----------------------------------------------------------------------------------------------------------
# Calculates the snapshot value of the first-derivative twist (a⁺, a⁻) in x, y and z directions.
function firstDerivativeTwist(ψ::Cuboid, a⁺::R, a⁻::R) where R<:Real

    f = ψ.syst.f
    L₁ = ψ.syst.L₁
    ν = ψ.syst.ν
    κ₅ = ψ.syst.κ₅

    # Define dynamic function to be used in lattice sum.
    function densityFunkX(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        firstDerivativeDensityX(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end
    function densityFunkY(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        firstDerivativeDensityY(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end
    function densityFunkZ(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂
        ϕᵣ₊₃ = nb.ϕᵣ₊₃

        firstDerivativeDensityZ(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, A⁰, z, ν, κ₅)
    end

    (latticeSiteNeighborSum(densityFunkX, cub, Float64), latticeSiteNeighborSum(densityFunkY, cub, Float64),
     latticeSiteNeighborSum(densityFunkZ, cub, Float64))
end

# -----------------------------------------------------------------------------------------------------------
# Calculate the snapshot value of the second-derivative (a⁺, a⁻) twist, in x, y and z directions.
function secondDerivativeTwist(ψ::Cuboid, a⁺::R, a⁻::R) where R<:Real

    f = ψ.syst.f
    L₁ = ψ.syst.L₁
    ν = ψ.syst.ν
    κ₅ = ψ.syst.κ₅

    # Define dynamic function to be used in lattice sum.
    function densityFunkX(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        secondDerivativeDensityX(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end
    function densityFunkY(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        secondDerivativeDensityY(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end
    function densityFunkZ(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂
        ϕᵣ₊₃ = nb.ϕᵣ₊₃

        secondDerivativeDensityZ(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, A⁰, z, ν, κ₅)
    end

    (latticeSiteNeighborSum(densityFunkX, cub, Float64), latticeSiteNeighborSum(densityFunkY, cub, Float64),
     latticeSiteNeighborSum(densityFunkZ, cub, Float64))
end

