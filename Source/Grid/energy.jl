############################################################################################################################
#                               Energy functions // SubCuboidModule
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# ---------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms. 
# This is a version of the energy above which is clearer to read, but should result in the same
# number but at the cost of some efficiency caused by calculating the square of u multiple times.
function fᵣ(ϕ::LatticeSite,nb::NearestNeighbors,c::SystConstants)
    energy = 0.0
	A₂ = ϕ.A[2]+two_pi*c.f*(ϕ.x-1)
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Complete kinetic term.
    Fₖ = 2*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1])
          + (ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)
    + c.κ₅*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3]))
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A[1])
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)
    + c.κ₅*((ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A[3])))
    Fₖt = -cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1]) - cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂) - c.κ₅*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3])
    # Potential energy term
    Fᵥ = ((ϕ.u⁺*ϕ.u⁻)^2*(2+c.ν*cos(2*(ϕ.θ⁻-ϕ.θ⁺))) 
          -(ϕ.u⁺)^2 + 0.5*(ϕ.u⁺)^4
          -(ϕ.u⁻)^2 + 0.5*(ϕ.u⁻)^4)
    # Anisotropi terms
    Fₐₙ = (c.ν+1)*(ϕ.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ - A₂) 
        + ϕ.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ - A₂) 
        - ϕ.u⁻*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - ϕ.A[1]) 
        - ϕ.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.A[1]))
    # Mixed gradient terms
    Fₘ = (c.ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*sin(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (ϕ.A[1]-A₂)) 
        - ϕᵣ₊₂.u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (ϕ.A[1] - A₂)) 
        + 2*ϕ.u⁺*ϕ.u⁻*sin(ϕ.θ⁻-ϕ.θ⁺) 
        + ϕ.u⁻*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ - A₂)    # Here there used to be a sign error
        - ϕ.u⁺*ϕᵣ₊₂.u⁻*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ - A₂)    # Here too
        + ϕ.u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - ϕ.A[1]) 
        - ϕ.u⁺*ϕᵣ₊₁.u⁻*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.A[1]))
    energy = Fₖt # + Fᵥ + Fₐₙ + Fₘ
end

# --------------------------------------------------------------------------------------------------
# Calculates the contribution to the maxwell term in the free energy from evaulating the sum at a
# specific site.
function maxwell(ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, g⁻²::Float64)
    return  ((ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2 + (ϕ.A[2] + ϕᵣ₊₂.A[3] - ϕᵣ₊₃.A[2] - ϕ.A[3])^2
                   + (ϕ.A[3] + ϕᵣ₊₃.A[1] - ϕᵣ₊₁.A[3] - ϕ.A[1])^2)*g⁻²
end


# --------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
# This is general ordinary lattice regularization energy without CP1.
function ΔE(ϕ′::LatticeSite, ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, s::SystConstants)
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
    A⁰ = two_pi*s.f*(h_pos-1)
    #const A⁰₊ = 2π*s.f*x
    #This section need to be re-checed vs E(ψ) by using test in Note 25 
    if h_pos == 1
        A⁰₋ = two_pi*s.f*(s.L₁-1)
    else
        A⁰₋ = two_pi*s.f*(h_pos-2)
    end
    
    # Normal kinetic terms
    δFₖ = 2*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺-ϕ′.A[1]) - ϕᵣ₋₁.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1]) - ϕᵣ₋₁.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]))
            + (ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺-ϕ′.A[2]-A⁰) - ϕᵣ₋₂.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-ϕ.A[2]-A⁰) - ϕᵣ₋₂.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰))
     + s.κ₅*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-ϕ′.A[3]) - ϕᵣ₋₃.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3])) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3]) - ϕᵣ₋₃.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3])))
           + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻-ϕ′.A[1]) - ϕᵣ₋₁.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1])) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A[1]) - ϕᵣ₋₁.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
            + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻-ϕ′.A[2]-A⁰) - ϕᵣ₋₂.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A[2]-A⁰)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-ϕ.A[2]-A⁰) - ϕᵣ₋₂.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A[2]-A⁰))
     + s.κ₅*((ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ′.θ⁻-ϕ′.A[3]) - ϕᵣ₋₃.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3])) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A[3]) - ϕᵣ₋₃.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3]))))
    # Rewrite of the above for the case when u⁺=1, u⁻=0
    δFₖt = -((cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺-ϕ′.A[1]) - cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1])) + (cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]))
           + (cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺-ϕ′.A[2]-A⁰) + cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰)) - (cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-ϕ.A[2]-A⁰) + cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰))
     + s.κ₅*((cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-ϕ′.A[3]) + cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3])) - (cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3]) + cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3]))))

    # Potential terms
    δFᵥ = (((ϕ′.u⁺*ϕ′.u⁻)^2*(2+s.ν*cos(2*(ϕ′.θ⁺-ϕ′.θ⁻))) - ϕ′.u⁺^2 + 0.5*ϕ′.u⁺^4 - ϕ′.u⁻^2 + 0.5*ϕ′.u⁻^4)
           - ((ϕ.u⁺*ϕ.u⁻)^2*(2+s.ν*cos(2*(ϕ.θ⁺-ϕ.θ⁻))) - ϕ.u⁺^2 + 0.5*ϕ.u⁺^4 - ϕ.u⁻^2 + 0.5*ϕ.u⁻^4))

    # Andreev-Bashkin terms
    δFₐₙ = (s.ν+1)*(ϕᵣ₊₂.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-(ϕ′.A[2]+A⁰)) - ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-(ϕ.A[2]+A⁰)))
        - ϕᵣ₊₁.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-ϕ′.A[1]) - ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₂.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-(ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-(ϕ.A[2]+A⁰))) 
        - ϕᵣ₊₁.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])))
    
    # Mixed gradient terms
    δFₘ = (s.ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*(sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ - (ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
            - sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ - (ϕ.A[1] - (ϕ.A[2]+A⁰)))) 
        +ϕᵣ₊₂.u⁺*(ϕ′.u⁻*sin(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-(ϕ′.A[2]+A⁰)) - ϕ.u⁻*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-(ϕ.A[2]+A⁰)))  # Old sign error
        +ϕᵣ₊₁.u⁺*(ϕ′.u⁻*sin(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-ϕ′.A[1]) - ϕ.u⁻*sin(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A[1])) 
        +ϕᵣ₋₂₊₁.u⁻*(ϕ′.u⁺*sin(ϕᵣ₋₂₊₁.θ⁻-ϕ′.θ⁺-(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰))) 
            - ϕ.u⁺*sin(ϕᵣ₋₂₊₁.θ⁻-ϕ.θ⁺-(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰)))) 
        +ϕᵣ₋₁₊₂.u⁺*(ϕ′.u⁻*sin(ϕ′.θ⁻-ϕᵣ₋₁₊₂.θ⁺-(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋))) 
            - ϕ.u⁻*sin(ϕ.θ⁻-ϕᵣ₋₁₊₂.θ⁺-(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋)))) 
        +ϕᵣ₋₂.u⁻*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)))  # Old sign error
        +ϕᵣ₋₁.u⁻*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1])) 
        -(ϕᵣ₊₂.u⁻*ϕᵣ₊₁.u⁺*(sin(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻ - (ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
            - sin(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻ - (ϕ.A[1] - (ϕ.A[2]+A⁰)))) 
        +ϕᵣ₊₂.u⁻*(ϕ′.u⁺*sin(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-(ϕ′.A[2]+A⁰)) - ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-(ϕ.A[2]+A⁰)))  # Old sign error
        +ϕᵣ₊₁.u⁻*(ϕ′.u⁺*sin(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-ϕ′.A[1]) - ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A[1])) 
        +ϕᵣ₋₂₊₁.u⁺*(ϕ′.u⁻*sin(ϕᵣ₋₂₊₁.θ⁺-ϕ′.θ⁻-(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰))) 
            - ϕ.u⁻*sin(ϕᵣ₋₂₊₁.θ⁺-ϕ.θ⁻-(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰)))) 
        +ϕᵣ₋₁₊₂.u⁻*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₁₊₂.θ⁻-(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋))) 
            - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₁₊₂.θ⁻-(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋)))) 
        +ϕᵣ₋₂.u⁺*(ϕ′.u⁻*sin(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*sin(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)))  # Old sign error
        +ϕᵣ₋₁.u⁺*(ϕ′.u⁻*sin(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁻*sin(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])))
        +2*(ϕ′.u⁺*ϕ′.u⁻*sin(ϕ′.θ⁻-ϕ′.θ⁺) - ϕ.u⁺*ϕ.u⁻*sin(ϕ.θ⁻-ϕ.θ⁺)))
    
    # Then calculate the Gauge field contribution
    ϕᵣ₊₁₋₂ = nnb.ϕᵣ₊₁₋₂
    ϕᵣ₊₁₋₃ = nnb.ϕᵣ₊₁₋₃
    ϕᵣ₋₁₊₃ = nnb.ϕᵣ₋₁₊₃
    ϕᵣ₊₂₋₃ = nnb.ϕᵣ₊₂₋₃
    ϕᵣ₋₂₊₃ = nnb.ϕᵣ₋₂₊₃
    # First contribution from current position
    #δE += ((ϕ′.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ′.A[2])^2 - (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])^2)*s.g⁻²
    δFₐ = maxwell(ϕ′, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, s.g⁻²) - maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, s.g⁻²)
    # Then from position r-x
    δFₐ += maxwell(ϕᵣ₋₁, ϕ′, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, s.g⁻²) - maxwell(ϕᵣ₋₁, ϕ, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, s.g⁻²)
    #δE += ((ϕᵣ₋₁.A[1] + ϕ′.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2 - (ϕᵣ₋₁.A[1] + ϕ.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2)*s.g⁻²
    # Then from position r-y
    δFₐ += maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ′, ϕᵣ₋₂₊₃, s.g⁻²) - maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ, ϕᵣ₋₂₊₃, s.g⁻²)
    #δE += ((ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ′.A[1] - ϕᵣ₋₂.A[2])^2 - (ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ.A[1] - ϕᵣ₋₂.A[2])^2)*s.g⁻²
    # Then we also have to include the position r-z
    δFₐ += maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ′, s.g⁻²) - maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ, s.g⁻²)

	δE = δFₐ + δFₖt# + δFᵥ + δFₐₙ + δFₘ 
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

