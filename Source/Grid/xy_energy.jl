############################################################################################################################
#                               Energy functions // SubCuboidModule
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# ---------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms. 
# This is a version of the energy above which is clearer to read, but should result in the same
# number but at the cost of some efficiency caused by calculating the square of u multiple times.
function fᵣ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
	A₂ = ϕ.A₂+two_pi*c.f*(ϕ.x-1)
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Complete kinetic term.
    Fₖ = 2*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A₁)
          + (ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)
    + c.κ₅*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A₃))
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A₁)
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)
    + c.κ₅*((ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A₃)))
    # Potential energy term
    Fᵥ = (0.5*(1+(1+c.ν)/2)*(ϕ.u⁺^4 + ϕ.u⁻^4) 
          + 0.5*(1-c.ν)*(ϕ.u⁺*ϕ.u⁻)^2*(2+cos(2*(ϕ.θ⁺-ϕ.θ⁻))) 
          - (ϕ.u⁺^2 + ϕ.u⁻^2))
    # Anisotropi terms
    Fₐₙ = (c.ν+1)*(ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - ϕ.A₁)
                 - ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - A₂)
                 + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - A₂)
                 - ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - ϕ.A₁))
    # Mixed gradient terms
    Fₘ = (1-c.ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-ϕ.A₁+A₂)
                -ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A₁)
                -ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-A₂)
                +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂+ϕ.A₁)
                -ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A₁)
                -ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-A₂)
                +2*ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻))
    energy = Fₖ + Fᵥ + Fₐₙ + Fₘ
end

# --------------------------------------------------------------------------------------------------
# Calculates the contribution to the maxwell term in the free energy from evaulating the sum at a
# specific site.
function maxwell(ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, g⁻²::Float64)
    return  ((ϕ.A₁ + ϕᵣ₊₁.A₂-ϕᵣ₊₂.A₁-ϕ.A₂)^2 + (ϕ.A₂ + ϕᵣ₊₂.A₃ - ϕᵣ₊₃.A₂ - ϕ.A₃)^2
                   + (ϕ.A₃ + ϕᵣ₊₃.A₁ - ϕᵣ₊₁.A₃ - ϕ.A₁)^2)*g⁻²
end


# --------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
# This is general ordinary lattice regularization energy without CP1.
function ΔE(ϕ′::LatticeSite, ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, c::SystConstants)
    δE::Float64 = 0.0
    h_pos = ϕ.x

	# Get neighbors
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃
	ϕᵣ₋₁ = nb.ϕᵣ₋₁
	ϕᵣ₋₂ = nb.ϕᵣ₋₂
    ϕᵣ₋₃ = nb.ϕᵣ₋₃

	ϕᵣ₋₁₊₂ = nnb.ϕᵣ₋₁₊₂
	ϕᵣ₋₂₊₁ = nnb.ϕᵣ₊₁₋₂
    
    # Calculate constant link variables
    A⁰ = two_pi*c.f*(h_pos-1)
    #const A⁰₊ = 2π*c.f*x
    #This section need to be re-checed vs E(ψ) by using test in Note 25 
    if h_pos == 1
        A⁰₋ = two_pi*c.f*(c.L₁-1)
    else
        A⁰₋ = two_pi*c.f*(h_pos-2)
    end
    
    # Normal kinetic terms
    δFₖ = 2*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺-ϕ′.A₁) - ϕᵣ₋₁.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A₁)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A₁) - ϕᵣ₋₁.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A₁))
            + (ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺-ϕ′.A₂-A⁰) - ϕᵣ₋₂.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A₂-A⁰)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-ϕ.A₂-A⁰) - ϕᵣ₋₂.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A₂-A⁰))
     + c.κ₅*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-ϕ′.A₃) - ϕᵣ₋₃.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A₃)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A₃) - ϕᵣ₋₃.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A₃)))
           + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻-ϕ′.A₁) - ϕᵣ₋₁.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A₁)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A₁) - ϕᵣ₋₁.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A₁))
            + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻-ϕ′.A₂-A⁰) - ϕᵣ₋₂.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A₂-A⁰)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-ϕ.A₂-A⁰) - ϕᵣ₋₂.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A₂-A⁰))
     + c.κ₅*((ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ′.θ⁻-ϕ′.A₃) - ϕᵣ₋₃.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A₃)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A₃) - ϕᵣ₋₃.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A₃))))

    # Potential terms
    δFᵥ = (0.5*(1+0.5*(1+c.ν))*((ϕ′.u⁺^4+ϕ′.u⁻^4) - (ϕ.u⁺^4+ϕ.u⁻^4)) 
         + 0.5*(1-c.ν)*((ϕ′.u⁺*ϕ′.u⁻)^2*(2+cos(2*(ϕ′.θ⁺-ϕ′.θ⁻))) - (ϕ.u⁺*ϕ.u⁻)^2*(2+cos(2*(ϕ.θ⁺-ϕ.θ⁻)))) 
         - (ϕ′.u⁺^2 + ϕ′.u⁻^2) + (ϕ.u⁺^2 + ϕ.u⁻^2))

    # Andreev-Bashkin terms
    δFₐₙ = (c.ν+1)*(ϕᵣ₊₂.u⁺*ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺ - (ϕ′.A₂+A⁰)) -  ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺ - (ϕ.A₂+A⁰)) 
                  + ϕ′.u⁺*ϕᵣ₋₂.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺ - (ϕᵣ₋₂.A₂+A⁰)) - ϕ.u⁺*ϕᵣ₋₂.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺ - (ϕᵣ₋₂.A₂+A⁰)) 
                  - ϕᵣ₊₁.u⁺*ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺ - ϕ′.A₁) + ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺ - ϕ.A₁)
                  - ϕ′.u⁺*ϕᵣ₋₁.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺ - ϕᵣ₋₁.A₁) + ϕ.u⁺*ϕᵣ₋₁.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺ - ϕᵣ₋₁.A₁)
                  + ϕᵣ₊₁.u⁻*ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻ - ϕ′.A₁) -  ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻ - ϕ.A₁) 
                  + ϕ′.u⁻*ϕᵣ₋₁.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻ - ϕᵣ₋₁.A₁) - ϕ.u⁻*ϕᵣ₋₁.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻ - ϕᵣ₋₁.A₁) 
                  - ϕᵣ₊₂.u⁻*ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻ - (ϕ′.A₂+A⁰)) + ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻ - (ϕ.A₂+A⁰))
                  - ϕ′.u⁻*ϕᵣ₋₂.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻ - (ϕᵣ₋₂.A₂+A⁰)) + ϕ.u⁻*ϕᵣ₋₂.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻ - (ϕᵣ₋₂.A₂+A⁰)))
    
    # Mixed gradient terms
    δFₘ = (1-c.ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-ϕ′.A₁+ϕ′.A₂+A⁰)
                 +ϕ′.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₁₊₂.θ⁻-ϕᵣ₋₁.A₁+ϕᵣ₋₁.A₂+A⁰₋)
                 +ϕᵣ₋₂₊₁.u⁺*ϕ′.u⁻*cos(ϕᵣ₋₂₊₁.θ⁺-ϕ′.θ⁻-ϕᵣ₋₂.A₁+ϕᵣ₋₂.A₂+A⁰)
                -(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-ϕ.A₁+ϕ.A₂+A⁰)
                 +ϕ.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₁₊₂.θ⁻-ϕᵣ₋₁.A₁+ϕᵣ₋₁.A₂+A⁰₋)
                 +ϕᵣ₋₂₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₋₂₊₁.θ⁺-ϕ.θ⁻-ϕᵣ₋₂.A₁+ϕᵣ₋₂.A₂+A⁰))
                 -ϕᵣ₊₁.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-ϕ′.A₁)-ϕ′.u⁺*ϕᵣ₋₁.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A₁)
                -(-ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A₁)-ϕ.u⁺*ϕᵣ₋₁.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A₁))
                 -ϕᵣ₊₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-ϕ′.A₂-A⁰)-ϕ′.u⁻*ϕᵣ₋₂.u⁺*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A₂-A⁰)
                -(-ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-ϕ.A₂-A⁰)-ϕ.u⁻*ϕᵣ₋₂.u⁺*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A₂-A⁰))
                 +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-ϕ′.A₂-A⁰+ϕ′.A₁)
                 +ϕ′.u⁺*ϕᵣ₋₂₊₁.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₂₊₁.θ⁻-ϕᵣ₋₂.A₂-A⁰+ϕᵣ₋₂.A₁)
                 +ϕᵣ₋₁₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺-ϕ′.θ⁻-ϕᵣ₋₁.A₂-A⁰₋+ϕᵣ₋₁.A₁)
                -(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-ϕ.A₂-A⁰+ϕ.A₁)
                 +ϕ.u⁺*ϕᵣ₋₂₊₁.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₂₊₁.θ⁻-ϕᵣ₋₂.A₂-A⁰+ϕᵣ₋₂.A₁)
                 +ϕᵣ₋₁₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺-ϕ.θ⁻-ϕᵣ₋₁.A₂-A⁰₋+ϕᵣ₋₁.A₁))
                 +2*ϕ′.u⁺*ϕ′.u⁻*cos(ϕ′.θ⁺-ϕ′.θ⁻)-2*ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻)
                 -ϕᵣ₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-ϕ′.A₂-A⁰)-ϕ′.u⁺*ϕᵣ₋₂.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A₂-A⁰)
                -(-ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-ϕ.A₂-A⁰)-ϕ.u⁺*ϕᵣ₋₂.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A₂-A⁰))
                 -ϕᵣ₊₁.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-ϕ′.A₁)-ϕ′.u⁻*ϕᵣ₋₁.u⁺*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A₁)
                -(-ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A₁)-ϕ.u⁻*ϕᵣ₋₁.u⁺*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A₁))) 
    
    # Then calculate the Gauge field contribution
    ϕᵣ₊₁₋₂ = nnb.ϕᵣ₊₁₋₂
    ϕᵣ₊₁₋₃ = nnb.ϕᵣ₊₁₋₃
    ϕᵣ₋₁₊₃ = nnb.ϕᵣ₋₁₊₃
    ϕᵣ₊₂₋₃ = nnb.ϕᵣ₊₂₋₃
    ϕᵣ₋₂₊₃ = nnb.ϕᵣ₋₂₊₃
    # First contribution from current position
    δFₐ = maxwell(ϕ′, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²) - maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²)
    # Then from position r-x
    δFₐ += maxwell(ϕᵣ₋₁, ϕ′, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²) - maxwell(ϕᵣ₋₁, ϕ, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²)
    # Then from position r-y
    δFₐ += maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ′, ϕᵣ₋₂₊₃, c.g⁻²) - maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ, ϕᵣ₋₂₊₃, c.g⁻²)
    # Then we also have to include the position r-z
    δFₐ += maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ′, c.g⁻²) - maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ, c.g⁻²)

	δE = δFₖ +  δFₐ + δFᵥ + δFₐₙ + δFₘ 
end

# Calculates the contribution to energy given by the sub-cuboid
function energy(sc::SubCuboid)
    g⁻² = sc.syst.g⁻²
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    E = 0.0
    
    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        E += fᵣ(ϕ,nb,sc.syst)
        E += maxwell(ϕ, nb.ϕᵣ₊₁, nb.ϕᵣ₊₂, nb.ϕᵣ₊₃, g⁻²)
    end
    E
end
# This function should be called from the worker process where the channel of the SubCuboid is located.
function energy(chan::RemoteChannel{Channel{SubCuboid}})
    sc = fetch(chan)
    energy(sc)
end

