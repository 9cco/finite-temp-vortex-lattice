############################################################################################################################
#                               Energy functions // SubCuboidModule
#__________________________________________________________________________________________________________________________#
############################################################################################################################
# The functions used to calculate energy and energy differences for doing local Markov-Chain Monte-Carlo updates. This 
# defines the energy of our model. In this case we have written the Γ₅ᵤ representation in the xy-basis such that 
# η⁺ → ηˣ and η⁻ → ηʸ.

# ---------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms. 
# This is a version of the energy above which is clearer to read, but should result in the same
# number but at the cost of some efficiency caused by calculating the square of u multiple times.
function fᵣ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
    ϕᵣ₊₁ = nb.ϕᵣ₊₁; ϕᵣ₋₁ = nb.ϕᵣ₋₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂; ϕᵣ₋₂ = nb.ϕᵣ₋₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃; ϕᵣ₋₃ = nb.ϕᵣ₋₃
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    Aᵣ₋₁_₁, Aᵣ₋₁_₂, Aᵣ₋₁_₃ = linkVariables(ϕᵣ₋₁, c)
    Aᵣ₋₂_₁, Aᵣ₋₂_₂, Aᵣ₋₂_₃ = linkVariables(ϕᵣ₋₂, c)

    # Complete kinetic term.
    Fₖ = 2*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-A₁)
          + (ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)
    + c.κ₅*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-A₃))
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-A₁)
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)
    + c.κ₅*((ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-A₃)))
    # Potential energy term
    Fᵥ = (0.5*(1+(1+c.ν)/2)*(ϕ.u⁺^4 + ϕ.u⁻^4) 
          + 0.5*(1-c.ν)*(ϕ.u⁺*ϕ.u⁻)^2*(2+cos(2*(ϕ.θ⁺-ϕ.θ⁻))) 
          - (ϕ.u⁺^2 + ϕ.u⁻^2))
    # Anisotropi terms
    Fₐₙ = (c.ν+1)*(ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - A₁)
                 - ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - A₂)
                 + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - A₂)
                 - ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - A₁))
    # Mixed gradient terms
    #Fₘ = c.κ*(1-c.ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-A₁+A₂)
    #            -ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-A₁)
    #            -ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-A₂)
    #            +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂+A₁)
    #            -ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-A₁)
    #            -ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-A₂)
    #            +2*ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻))
    # Symmetrized mixed gradient terms
    Fₘ = (ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-A₁+A₂)
                +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂+A₁))
    Fₘ_C₄ = -(ϕᵣ₋₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₋₁.θ⁺-ϕᵣ₊₂.θ⁻ + A₂ + Aᵣ₋₁_₁)
                +ϕᵣ₋₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₋₁.θ⁻-ϕᵣ₊₂.θ⁺ + A₂ + Aᵣ₋₁_₁))
    Fₘ_C₄² = (ϕᵣ₋₁.u⁺*ϕᵣ₋₂.u⁻*cos(ϕᵣ₋₁.θ⁺-ϕᵣ₋₂.θ⁻ - (Aᵣ₋₂_₂ - Aᵣ₋₁_₁))
        + ϕᵣ₋₁.u⁻*ϕᵣ₋₂.u⁺*cos(ϕᵣ₋₁.θ⁻-ϕᵣ₋₂.θ⁺ - (Aᵣ₋₂_₂ - Aᵣ₋₁_₁)))
    Fₘ_C₄³ = -(ϕᵣ₊₁.u⁺*ϕᵣ₋₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₋₂.θ⁻ - (Aᵣ₋₂_₂ + A₁))
        + ϕᵣ₊₁.u⁻*ϕᵣ₋₂.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕᵣ₋₂.θ⁺ - (Aᵣ₋₂_₂ + A₁)))
    Fₘ_symm = c.κ*(1-c.ν)*(Fₘ + Fₘ_C₄ + Fₘ_C₄² + Fₘ_C₄³)/4

    energy = Fₖ + Fᵥ + Fₐₙ + Fₘ_symm
end

# --------------------------------------------------------------------------------------------------
# Calculates the contribution to the maxwell term in the free energy from evaulating the sum at a
# specific site. This includes both the constant as well as the fluctuating gauge field.
function maxwell(ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, s::SystConstants)
    A₁, A₂, A₃ = ϕ.A₁, ϕ.A₂, ϕ.A₃ #linkVariables(ϕ, s)
    Aʳ⁺¹₂ = ϕᵣ₊₁.A₂ #linkVariableY(ϕᵣ₊₁, s)
    Aʳ⁺¹₃ = ϕᵣ₊₁.A₃ #linkVariableZ(ϕᵣ₊₁, s)
    Aʳ⁺²₁ = ϕᵣ₊₂.A₁ #linkVariableX(ϕᵣ₊₂, s)
    Aʳ⁺²₃ = ϕᵣ₊₂.A₃ #linkVariableZ(ϕᵣ₊₂, s)
    Aʳ⁺³₁ = ϕᵣ₊₃.A₁ #linkVariableX(ϕᵣ₊₃, s)
    Aʳ⁺³₂ = ϕᵣ₊₃.A₂ #linkVariableY(ϕᵣ₊₃, s)
    return  ((A₁ + Aʳ⁺¹₂- Aʳ⁺²₁ -A₂)^2 + (A₂ + Aʳ⁺²₃ - Aʳ⁺³₂ - A₃)^2
                   + (A₃ + Aʳ⁺³₁ - Aʳ⁺¹₃ - A₁)^2)*s.g⁻²
end


# --------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
# This is general ordinary lattice regularization energy without CP1.
function ΔE(ϕ′::LatticeSite, ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, c::SystConstants)
    δE::Float64 = 0.0

	# Get neighbors
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃
	ϕᵣ₋₁ = nb.ϕᵣ₋₁
	ϕᵣ₋₂ = nb.ϕᵣ₋₂
    ϕᵣ₋₃ = nb.ϕᵣ₋₃

	ϕᵣ₋₁₊₂ = nnb.ϕᵣ₋₁₊₂
    ϕᵣ₋₁₋₂ = nnb.ϕᵣ₋₁₋₂
	ϕᵣ₊₁₋₂ = nnb.ϕᵣ₊₁₋₂
    ϕᵣ₊₁₊₂ = nnb.ϕᵣ₊₁₊₂
    ϕᵣ₊₁₋₃ = nnb.ϕᵣ₊₁₋₃
    ϕᵣ₋₁₊₃ = nnb.ϕᵣ₋₁₊₃
    ϕᵣ₊₂₋₃ = nnb.ϕᵣ₊₂₋₃
    ϕᵣ₋₂₊₃ = nnb.ϕᵣ₋₂₊₃
    
    # Calculate constant link variables
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    A₁′, A₂′, A₃′ = linkVariables(ϕ′, c)
    Aᵣ₋₁_₁ = linkVariableX(ϕᵣ₋₁, c)
    Aᵣ₋₁_₂ = linkVariableY(ϕᵣ₋₁, c)
    Aᵣ₋₂_₁ = linkVariableX(ϕᵣ₋₂, c)
    Aᵣ₋₂_₂ = linkVariableY(ϕᵣ₋₂, c)
    Aᵣ₋₃_₃ = linkVariableZ(ϕᵣ₋₃, c)
    Aᵣ₋₁₋₂_₁ = linkVariableX(ϕᵣ₋₁₋₂, c)
    Aᵣ₋₁₋₂_₂ = linkVariableY(ϕᵣ₋₁₋₂, c)

    Aᵣ₊₁_₂ = linkVariableY(ϕᵣ₊₁, c)
    Aᵣ₊₂_₁ = linkVariableX(ϕᵣ₊₂, c)
    Aᵣ₊₁₋₂_₂ = linkVariableY(ϕᵣ₊₁₋₂, c)
    Aᵣ₋₁₊₂_₁ = linkVariableX(ϕᵣ₋₁₊₂, c)

    
    # Normal kinetic terms
    δFₖ = 2*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺-A₁′) - ϕᵣ₋₁.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-Aᵣ₋₁_₁)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-A₁) - ϕᵣ₋₁.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-Aᵣ₋₁_₁))
            + (ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺-A₂′) - ϕᵣ₋₂.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-Aᵣ₋₂_₂)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂) - ϕᵣ₋₂.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-Aᵣ₋₂_₂))
     + c.κ₅*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-A₃′) - ϕᵣ₋₃.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-Aᵣ₋₃_₃)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-A₃) - ϕᵣ₋₃.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-Aᵣ₋₃_₃)))
           + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻-A₁′) - ϕᵣ₋₁.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-Aᵣ₋₁_₁)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-A₁) - ϕᵣ₋₁.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-Aᵣ₋₁_₁))
            + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻-A₂′) - ϕᵣ₋₂.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-Aᵣ₋₂_₂)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂) - ϕᵣ₋₂.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-Aᵣ₋₂_₂))
     + c.κ₅*((ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ′.θ⁻-A₃′) - ϕᵣ₋₃.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₃.θ⁻-Aᵣ₋₃_₃)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-A₃) - ϕᵣ₋₃.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₃.θ⁻-Aᵣ₋₃_₃))))

    # Potential terms
    δFᵥ = (0.5*(1+0.5*(1+c.ν))*((ϕ′.u⁺^4+ϕ′.u⁻^4) - (ϕ.u⁺^4+ϕ.u⁻^4)) 
         + 0.5*(1-c.ν)*((ϕ′.u⁺*ϕ′.u⁻)^2*(2+cos(2*(ϕ′.θ⁺-ϕ′.θ⁻))) - (ϕ.u⁺*ϕ.u⁻)^2*(2+cos(2*(ϕ.θ⁺-ϕ.θ⁻)))) 
         - (ϕ′.u⁺^2 + ϕ′.u⁻^2) + (ϕ.u⁺^2 + ϕ.u⁻^2))

    # Andreev-Bashkin terms
    δFₐₙ = (c.ν+1)*(ϕᵣ₊₂.u⁺*ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺ - (A₂′)) -  ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺ - (A₂)) 
                  + ϕ′.u⁺*ϕᵣ₋₂.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺ - (Aᵣ₋₂_₂)) - ϕ.u⁺*ϕᵣ₋₂.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺ - (Aᵣ₋₂_₂)) 
                  - ϕᵣ₊₁.u⁺*ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺ - A₁′) + ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺ - A₁)
                  - ϕ′.u⁺*ϕᵣ₋₁.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺ - Aᵣ₋₁_₁) + ϕ.u⁺*ϕᵣ₋₁.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺ - Aᵣ₋₁_₁)
                  + ϕᵣ₊₁.u⁻*ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻ - A₁′) -  ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻ - A₁) 
                  + ϕ′.u⁻*ϕᵣ₋₁.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻ - Aᵣ₋₁_₁) - ϕ.u⁻*ϕᵣ₋₁.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻ - Aᵣ₋₁_₁) 
                  - ϕᵣ₊₂.u⁻*ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻ - (A₂′)) + ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻ - (A₂))
                  - ϕ′.u⁻*ϕᵣ₋₂.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻ - (Aᵣ₋₂_₂)) + ϕ.u⁻*ϕᵣ₋₂.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻ - (Aᵣ₋₂_₂)))
    
    # Mixed gradient terms
#    δFₘ = c.κ*(1-c.ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-A₁′+A₂′)
#                 +ϕ′.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₁₊₂.θ⁻-Aᵣ₋₁_₁+Aᵣ₋₁_₂)
#                 +ϕᵣ₊₁₋₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₁₋₂.θ⁺-ϕ′.θ⁻-Aᵣ₋₂_₁+Aᵣ₋₂_₂)
#                -(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-A₁+A₂)
#                 +ϕ.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₁₊₂.θ⁻-Aᵣ₋₁_₁+Aᵣ₋₁_₂)
#                 +ϕᵣ₊₁₋₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁₋₂.θ⁺-ϕ.θ⁻-Aᵣ₋₂_₁+Aᵣ₋₂_₂))
#                 -ϕᵣ₊₁.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-A₁′)-ϕ′.u⁺*ϕᵣ₋₁.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-Aᵣ₋₁_₁)
#                -(-ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-A₁)-ϕ.u⁺*ϕᵣ₋₁.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-Aᵣ₋₁_₁))
#                 -ϕᵣ₊₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-A₂′)-ϕ′.u⁻*ϕᵣ₋₂.u⁺*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-Aᵣ₋₂_₂)
#                -(-ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-A₂)-ϕ.u⁻*ϕᵣ₋₂.u⁺*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-Aᵣ₋₂_₂))
#                 +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂′+A₁′)
#                 +ϕ′.u⁺*ϕᵣ₊₁₋₂.u⁻*cos(ϕ′.θ⁺-ϕᵣ₊₁₋₂.θ⁻-Aᵣ₋₂_₂+Aᵣ₋₂_₁)
#                 +ϕᵣ₋₁₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺-ϕ′.θ⁻-Aᵣ₋₁_₂+Aᵣ₋₁_₁)
#                -(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂+A₁)
#                 +ϕ.u⁺*ϕᵣ₊₁₋₂.u⁻*cos(ϕ.θ⁺-ϕᵣ₊₁₋₂.θ⁻-Aᵣ₋₂_₂+Aᵣ₋₂_₁)
#                 +ϕᵣ₋₁₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺-ϕ.θ⁻-Aᵣ₋₁_₂+Aᵣ₋₁_₁))
#                 +2*ϕ′.u⁺*ϕ′.u⁻*cos(ϕ′.θ⁺-ϕ′.θ⁻)-2*ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻)
#                 -ϕᵣ₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-A₂′)-ϕ′.u⁺*ϕᵣ₋₂.u⁻*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-Aᵣ₋₂_₂)
#                -(-ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-A₂)-ϕ.u⁺*ϕᵣ₋₂.u⁻*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-Aᵣ₋₂_₂))
#                 -ϕᵣ₊₁.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-A₁′)-ϕ′.u⁻*ϕᵣ₋₁.u⁺*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-Aᵣ₋₁_₁)
#                -(-ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-A₁)-ϕ.u⁻*ϕᵣ₋₁.u⁺*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-Aᵣ₋₁_₁))) 

    # Symmetrized mixed gradient terms
    δFₘ_sy = c.κ*(1-c.ν)/4*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (A₁′ - A₂′))
                      - ϕᵣ₋₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₋₁.θ⁺ - ϕᵣ₊₂.θ⁻ + (A₂′ + Aᵣ₋₁_₁))
                      - ϕᵣ₊₁.u⁺*ϕᵣ₋₂.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕᵣ₋₂.θ⁻ + (-Aᵣ₋₂_₂ - A₁′))

                      + ϕ′.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ′.θ⁺ - ϕᵣ₋₁₊₂.θ⁻ - (Aᵣ₋₁_₁ - Aᵣ₋₁_₂))
                      - ϕ′.u⁺*ϕᵣ₋₁₋₂.u⁻*cos(ϕ′.θ⁺ - ϕᵣ₋₁₋₂.θ⁻ + (-Aᵣ₋₁₋₂_₂ - Aᵣ₋₁_₁))
                      
                      + ϕᵣ₊₁₋₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₁₋₂.θ⁺ - ϕ′.θ⁻ - (Aᵣ₋₂_₁ - Aᵣ₋₂_₂))
                      - ϕᵣ₋₁₋₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₋₁₋₂.θ⁺ - ϕ′.θ⁻ + (Aᵣ₋₂_₂ + Aᵣ₋₁₋₂_₁))

                      - ϕ′.u⁺*ϕᵣ₊₁₊₂.u⁻*cos(ϕ′.θ⁺ - ϕᵣ₊₁₊₂.θ⁻ + (Aᵣ₊₁_₂ + A₁′))
                      + ϕ′.u⁺*ϕᵣ₊₁₋₂.u⁻*cos(ϕ′.θ⁺ - ϕᵣ₊₁₋₂.θ⁻ - (Aᵣ₊₁₋₂_₂ - A₁′))

                      + ϕᵣ₋₁₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺ - ϕ′.θ⁻ - (A₂′ - Aᵣ₋₁₊₂_₁))
                      - ϕᵣ₊₁₊₂.u⁺*ϕ′.u⁻*cos(ϕᵣ₊₁₊₂.θ⁺ - ϕ′.θ⁻ + (-A₂′ - Aᵣ₊₂_₁))


                      -(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (A₁ - A₂))
                      - ϕᵣ₋₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₋₁.θ⁺ - ϕᵣ₊₂.θ⁻ + (A₂ + Aᵣ₋₁_₁))
                      - ϕᵣ₊₁.u⁺*ϕᵣ₋₂.u⁻*cos(ϕᵣ₊₁.θ⁺ - ϕᵣ₋₂.θ⁻ + (-Aᵣ₋₂_₂ - A₁))

                      + ϕ.u⁺*ϕᵣ₋₁₊₂.u⁻*cos(ϕ.θ⁺ - ϕᵣ₋₁₊₂.θ⁻ - (Aᵣ₋₁_₁ - Aᵣ₋₁_₂))
                      - ϕ.u⁺*ϕᵣ₋₁₋₂.u⁻*cos(ϕ.θ⁺ - ϕᵣ₋₁₋₂.θ⁻ + (-Aᵣ₋₁₋₂_₂ - Aᵣ₋₁_₁))
                      
                      + ϕᵣ₊₁₋₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁₋₂.θ⁺ - ϕ.θ⁻ - (Aᵣ₋₂_₁ - Aᵣ₋₂_₂))
                      - ϕᵣ₋₁₋₂.u⁺*ϕ.u⁻*cos(ϕᵣ₋₁₋₂.θ⁺ - ϕ.θ⁻ + (Aᵣ₋₂_₂ + Aᵣ₋₁₋₂_₁))

                      - ϕ.u⁺*ϕᵣ₊₁₊₂.u⁻*cos(ϕ.θ⁺ - ϕᵣ₊₁₊₂.θ⁻ + (Aᵣ₊₁_₂ + A₁))
                      + ϕ.u⁺*ϕᵣ₊₁₋₂.u⁻*cos(ϕ.θ⁺ - ϕᵣ₊₁₋₂.θ⁻ - (Aᵣ₊₁₋₂_₂ - A₁))

                      + ϕᵣ₋₁₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₋₁₊₂.θ⁺ - ϕ.θ⁻ - (A₂ - Aᵣ₋₁₊₂_₁))
                      - ϕᵣ₊₁₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁₊₂.θ⁺ - ϕ.θ⁻ + (-A₂ - Aᵣ₊₂_₁)))


                      + ϕᵣ₊₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (A₁′ - A₂′))
                      - ϕᵣ₋₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₋₁.θ⁻ - ϕᵣ₊₂.θ⁺ + (A₂′ + Aᵣ₋₁_₁))
                      - ϕᵣ₊₁.u⁻*ϕᵣ₋₂.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕᵣ₋₂.θ⁺ + (-Aᵣ₋₂_₂ - A₁′))

                      + ϕ′.u⁻*ϕᵣ₋₁₊₂.u⁺*cos(ϕ′.θ⁻ - ϕᵣ₋₁₊₂.θ⁺ - (Aᵣ₋₁_₁ - Aᵣ₋₁_₂))
                      - ϕ′.u⁻*ϕᵣ₋₁₋₂.u⁺*cos(ϕ′.θ⁻ - ϕᵣ₋₁₋₂.θ⁺ + (-Aᵣ₋₁₋₂_₂ - Aᵣ₋₁_₁))
                      
                      + ϕᵣ₊₁₋₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₁₋₂.θ⁻ - ϕ′.θ⁺ - (Aᵣ₋₂_₁ - Aᵣ₋₂_₂))
                      - ϕᵣ₋₁₋₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₋₁₋₂.θ⁻ - ϕ′.θ⁺ + (Aᵣ₋₂_₂ + Aᵣ₋₁₋₂_₁))

                      - ϕ′.u⁻*ϕᵣ₊₁₊₂.u⁺*cos(ϕ′.θ⁻ - ϕᵣ₊₁₊₂.θ⁺ + (Aᵣ₊₁_₂ + A₁′))
                      + ϕ′.u⁻*ϕᵣ₊₁₋₂.u⁺*cos(ϕ′.θ⁻ - ϕᵣ₊₁₋₂.θ⁺ - (Aᵣ₊₁₋₂_₂ - A₁′))

                      + ϕᵣ₋₁₊₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₋₁₊₂.θ⁻ - ϕ′.θ⁺ - (A₂′ - Aᵣ₋₁₊₂_₁))
                      - ϕᵣ₊₁₊₂.u⁻*ϕ′.u⁺*cos(ϕᵣ₊₁₊₂.θ⁻ - ϕ′.θ⁺ + (-A₂′ - Aᵣ₊₂_₁))


                      -(ϕᵣ₊₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (A₁ - A₂))
                      - ϕᵣ₋₁.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₋₁.θ⁻ - ϕᵣ₊₂.θ⁺ + (A₂ + Aᵣ₋₁_₁))
                      - ϕᵣ₊₁.u⁻*ϕᵣ₋₂.u⁺*cos(ϕᵣ₊₁.θ⁻ - ϕᵣ₋₂.θ⁺ + (-Aᵣ₋₂_₂ - A₁))

                      + ϕ.u⁻*ϕᵣ₋₁₊₂.u⁺*cos(ϕ.θ⁻ - ϕᵣ₋₁₊₂.θ⁺ - (Aᵣ₋₁_₁ - Aᵣ₋₁_₂))
                      - ϕ.u⁻*ϕᵣ₋₁₋₂.u⁺*cos(ϕ.θ⁻ - ϕᵣ₋₁₋₂.θ⁺ + (-Aᵣ₋₁₋₂_₂ - Aᵣ₋₁_₁))
                      
                      + ϕᵣ₊₁₋₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁₋₂.θ⁻ - ϕ.θ⁺ - (Aᵣ₋₂_₁ - Aᵣ₋₂_₂))
                      - ϕᵣ₋₁₋₂.u⁻*ϕ.u⁺*cos(ϕᵣ₋₁₋₂.θ⁻ - ϕ.θ⁺ + (Aᵣ₋₂_₂ + Aᵣ₋₁₋₂_₁))

                      - ϕ.u⁻*ϕᵣ₊₁₊₂.u⁺*cos(ϕ.θ⁻ - ϕᵣ₊₁₊₂.θ⁺ + (Aᵣ₊₁_₂ + A₁))
                      + ϕ.u⁻*ϕᵣ₊₁₋₂.u⁺*cos(ϕ.θ⁻ - ϕᵣ₊₁₋₂.θ⁺ - (Aᵣ₊₁₋₂_₂ - A₁))

                      + ϕᵣ₋₁₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₋₁₊₂.θ⁻ - ϕ.θ⁺ - (A₂ - Aᵣ₋₁₊₂_₁))
                      - ϕᵣ₊₁₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁₊₂.θ⁻ - ϕ.θ⁺ + (-A₂ - Aᵣ₊₂_₁))))
    
    # Then calculate the Gauge field contribution⁺⁻
    # First contribution from current position
    δFₐ = maxwell(ϕ′, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c) - maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c)
    # Then from position r-x
    δFₐ += maxwell(ϕᵣ₋₁, ϕ′, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c) - maxwell(ϕᵣ₋₁, ϕ, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c)
    # Then from position r-y
    δFₐ += maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ′, ϕᵣ₋₂₊₃, c) - maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ, ϕᵣ₋₂₊₃, c)
    # Then we also have to include the position r-z
    δFₐ += maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ′, c) - maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ, c)

	δE = δFₖ +  δFₐ + δFᵥ + δFₐₙ + δFₘ_sy 
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
        E += maxwell(ϕ, nb.ϕᵣ₊₁, nb.ϕᵣ₊₂, nb.ϕᵣ₊₃, sc.syst)
    end
    E
end
# This function should be called from the worker process where the channel of the SubCuboid is located.
function energy(chan::RemoteChannel{Channel{SubCuboid}})
    sc = fetch(chan)
    energy(sc)
end

