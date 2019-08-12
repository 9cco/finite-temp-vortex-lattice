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
# Calculates density of the first derivative of the free energy w.r.t. a twist in x and y direction.
# Since there are two phases the function calculates for a general twist coverned by a⁺ and a⁻
function firstDerivativeDensity(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    fk_x = 2*(ϕᵣ₊₁.u⁺*ϕ.u⁺*sin(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-A₁)a⁺ 
          + ϕᵣ₊₁.u⁻*ϕ.u⁻*sin(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-A₁)a⁻)
    fmc_x = (ν+1)*(ϕᵣ₊₁.u⁺*ϕ.u⁻*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))
               + ϕᵣ₊₁.u⁻*ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁)*(a⁻-x*(a⁺-a⁻))
               - x*(a⁺-a⁻)*ϕᵣ₊₂.u⁺*ϕ.u⁻*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)
               + x*(a⁺-a⁻)*ϕᵣ₊₂.u⁻*ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂))
    fmg_x = (1-ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (A₁ - A₂))*(a⁺+x*(a⁺-a⁻)) 
                 + ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻)) 
                 + x*(a⁺-a⁻)*(ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂) + ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻))
                 + ϕᵣ₊₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (A₁ - A₂))*(-a⁻+x*(a⁺-a⁻)) 
                 + ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻)) 
                 + x*(a⁺-a⁻)*(ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂) + ϕ.u⁻*ϕ.u⁺*cos(ϕ.θ⁻-ϕ.θ⁺)))
    fv_x = -ν*(ϕ.u⁺*ϕ.u⁻)^2*x*(a⁺-a⁻)*sin(2(ϕ.θ⁺ - ϕ.θ⁻))


    fk_y = 2*(ϕᵣ₊₂.u⁺*ϕ.u⁺*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)*a⁺ 
              + ϕᵣ₊₂.u⁻*ϕ.u⁻*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)a⁻)
    fmc_y = -(ν+1)*(ϕᵣ₊₂.u⁺*ϕ.u⁻*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))
                    + ϕᵣ₊₂.u⁻*ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂)*(a⁻-y*(a⁺-a⁻))
               - y*(a⁺-a⁻)*ϕᵣ₊₁.u⁺*ϕ.u⁻*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)
               + y*(a⁺-a⁻)*ϕᵣ₊₁.u⁻*ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁))
    fmg_y = (1-ν)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺ - ϕᵣ₊₁.θ⁻ - (A₂ - A₁))*(-a⁺+y*(a⁺-a⁻)) 
                   + ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻)) 
                   + y*(a⁺-a⁻)*(ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁) + ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻))
                   + ϕᵣ₊₂.u⁻*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₂.θ⁻ - ϕᵣ₊₁.θ⁺ - (A₂ - A₁))*(a⁻+y*(a⁺-a⁻)) 
                   + ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻)) 
                   + y*(a⁺-a⁻)*(ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁) + ϕ.u⁻*ϕ.u⁺*cos(ϕ.θ⁻-ϕ.θ⁺)))
    fv_y = -ν*(ϕ.u⁺*ϕ.u⁻)^2*y*(a⁺-a⁻)*sin(2(ϕ.θ⁺ - ϕ.θ⁻))

    fk_x+fmc_x+fmg_x+fv_x, fk_y+fmc_y+fmg_y+fv_y
end

# -----------------------------------------------------------------------------------------------------------
# Calculates the density of second derivative of the free energy w.r.t. a twist in x and y direction
# for a general (a⁺, a⁻) twist.
function secondDerivativeDensity(a⁺::R, a⁻::R, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, A::P, x::I, y::I, ν::P, κ₅::P) 
    where {P<:Real, I<:Int, R<:Real}

    A₂ = ϕ.A₂ + A
    A₁ = ϕ.A₁

    ddk_x = 2*(ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-A₁)(a⁺)^2
               + ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-A₁)(a⁻)^2)
    ddmc_x = (ν+1)*(ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))^2
               + ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻))^2
               - x^2*(a⁺-a⁻)^2*ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)
               - x^2*(a⁺-a⁻)^2*ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂))
    ddmg_x = (ν-1)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*sin(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (A₁ - A₂))*(a⁺+x*(a⁺-a⁻))^2 
                 + ϕᵣ₊₁.u⁺*ϕ.u⁻*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)*(a⁺+x*(a⁺-a⁻))^2 
                 + x^2*(a⁺-a⁻)^2*(ϕᵣ₊₂.u⁺*ϕ.u⁻*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂) + ϕ.u⁺*ϕ.u⁻*sin(ϕ.θ⁺-ϕ.θ⁻))
                 - ϕᵣ₊₁.u⁻*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (A₁ - A₂))*(-a⁻+x*(a⁺-a⁻))^2 
                 - ϕᵣ₊₁.u⁻*ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁)*(-a⁻+x*(a⁺-a⁻))^2 
                 - x^2*(a⁺-a⁻)^2*(ϕᵣ₊₂.u⁻*ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂) + ϕ.u⁻*ϕ.u⁺*sin(ϕ.θ⁻-ϕ.θ⁺)))
    ddv_x = -ν*(ϕ.u⁺*ϕ.u⁻)^2*x^2*(a⁺-a⁻)^2*cos(2(ϕ.θ⁺ - ϕ.θ⁻))


    ddk_y = 2*(ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)(a⁺)^2
               + ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)(a⁻)^2)
    ddmc_y = -(ν+1)*(ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))^2
                     + ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻))^2
               - y^2*(a⁺-a⁻)^2*ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁)
               - y^2*(a⁺-a⁻)^2*ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁))
    ddmg_y = (ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*sin(ϕᵣ₊₂.θ⁺ - ϕᵣ₊₁.θ⁻ - (A₂ - A₁))*(a⁺-y*(a⁺-a⁻))^2 
                    + ϕᵣ₊₂.u⁺*ϕ.u⁻*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁻ - A₂)*(a⁺+y*(a⁺-a⁻))^2 
                 + y^2*(a⁺-a⁻)^2*(ϕᵣ₊₁.u⁺*ϕ.u⁻*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - A₁) + ϕ.u⁺*ϕ.u⁻*sin(ϕ.θ⁺-ϕ.θ⁻))
                 - ϕᵣ₊₂.u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₂.θ⁻ - ϕᵣ₊₁.θ⁺ - (A₂ - A₁))*(-a⁻-y*(a⁺-a⁻))^2 
                 - ϕᵣ₊₂.u⁻*ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - A₂)*(-a⁻+y*(a⁺-a⁻))^2 
                 - y^2*(a⁺-a⁻)^2*(ϕᵣ₊₁.u⁻*ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - A₁) + ϕ.u⁻*ϕ.u⁺*sin(ϕ.θ⁻-ϕ.θ⁺)))
    ddv_y = -ν*(ϕ.u⁺*ϕ.u⁻)^2*y^2*(a⁺-a⁻)^2*cos(2(ϕ.θ⁺ - ϕ.θ⁻))

    ddk_x+ddmc_x+ddmg_x+ddv_x, ddk_y+dddmc_y+ddmg_y+ddv_y
end

# -----------------------------------------------------------------------------------------------------------
# Calculates the snapshot value of the first-derivative twist (a⁺, a⁻) in x and y directions.
function firstDerivativeTwist(ψ::Cuboid, a⁺::R, a⁻::R) where R<:Real

    f = ψ.syst.f
    L₁ = ψ.syst.L₁
    ν = ψ.syst.ν
    κ₅ = ψ.syst.κ₅

    # Define dynamic function to be used in lattice sum.
    function densityFunk(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        firstDerivativeDensity(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end

    (dH_x, dH_y) = latticeSiteNeighborSum(densityFunk, cub, Float64)
end

# -----------------------------------------------------------------------------------------------------------
# Calculate the snapshot value of the second-derivative (a⁺, a⁻) twist, in x and y directions.
function secondDerivativeTwist(ψ::Cuboid, a⁺::R, a⁻::R) where R<:Real

    f = ψ.syst.f
    L₁ = ψ.syst.L₁
    ν = ψ.syst.ν
    κ₅ = ψ.syst.κ₅

    # Define dynamic function to be used in lattice sum.
    function densityFunk(ϕ::LatticeSite, nb::NearestNeighbors, x::I, y::I, z::I) where I<:Int

        A⁰ = two_pi*f*(x-1)
        ϕᵣ₊₁ = nb.ϕᵣ₊₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂

        secondDerivativeDensity(a⁺, a⁻, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, A⁰, x, y, ν, κ₅)
    end

    (d²H_x, d²H_y) = latticeSiteNeighborSum(densityFunk, cub, Float64)
end

