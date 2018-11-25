####################################################################################################
#                            Energy functions
#
####################################################################################################

# ---------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms. 
# This is a version of the energy above which is clearer to read, but should result in the same
# number but at the cost of some efficiency caused by calculating the square of u multiple times.
function fᵣ(ϕ::LatticeSite,nb::NearestNeighbors,h_pos::Int64,c::SystConstants)
    energy = 0.0
	A₂ = ϕ.A[2]+two_pi*c.f*(h_pos-1)
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Complete kinetic term.
    Fₖ = 2*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1])
          + (ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)
     + c.κ₅*(ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3])
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A[1])
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)
     + c.κ₅*(ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A[3]))
    # Potential energy term
    Fᵥ = ((ϕ.u⁺*ϕ.u⁻)^2*(2+c.ν*cos(2*(ϕ.θ⁻-ϕ.θ⁺))) 
          -(1-1/c.β)*(ϕ.u⁺)^2 + 0.5*(ϕ.u⁺)^4
          -(1-1/c.β)*(ϕ.u⁻)^2 + 0.5*(ϕ.u⁻)^4)
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
    energy = Fₖ + Fᵥ #+ Fₐₙ + Fₘ
end

# ---------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms. This gives
# the CP1 constrained theory.
function fᵣCP1(ϕ::LatticeSite,nb::NearestNeighbors,h_pos::Int64,c::SystConstants)
    energy = 0.0
	A₂ = ϕ.A[2]+two_pi*c.f*(h_pos-1)
	ϕᵣ₊₁ = nb.ϕᵣ₊₁
	ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Kinetic energy Fₖ
    Fₖ = -2*(ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺ - ϕ.A[1]) 
        + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺ - A₂) 
        + ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻ - ϕ.A[1]) 
        + ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻ - A₂))
    # z-direction term
    Fₖ += 2*c.κ₅*(ϕ.u⁺^2 + ϕ.u⁻^2
        - ϕᵣ₊₃.u⁺*ϕ.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺ - ϕ.A[3])
        - ϕᵣ₊₃.u⁻*ϕ.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻ - ϕ.A[3]))
    # Potential energy Fᵥ
    Fᵥ = ϕ.u⁺^2*ϕ.u⁻^2*(1+c.ν*cos(2*(ϕ.θ⁺-ϕ.θ⁻)))
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
    energy = Fₖ + Fᵥ + Fₐₙ + Fₘ
end

function maxwell(ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₃::LatticeSite, g⁻²::Float64)
    return  ((ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2 + (ϕ.A[2] + ϕᵣ₊₂.A[3] - ϕᵣ₊₃.A[2] - ϕ.A[3])^2
                   + (ϕ.A[3] + ϕᵣ₊₃.A[1] - ϕᵣ₊₁.A[3] - ϕ.A[1])^2)*g⁻²
end


# --------------------------------------------------------------------------------------------------
# Loops over all positions of the lattice of a state and calculates the total energy from the
# Higgs-field terms using the function fᵣ() + the energy from the gauge field.
function E(ψ::State)
    energy = 0.0
    L = ψ.consts.L
	g⁻² = ψ.consts.g⁻²
    
    # Contribution from the bulk of lattice sites
    for z_pos=1:ψ.consts.L₃, h_pos=1:L, v_pos=1:L
        ϕ = ψ.lattice[v_pos,h_pos,z_pos]          # Lattice site at position r
		ϕᵣ₊₁ = ψ.nb[v_pos,h_pos,z_pos].ϕᵣ₊₁
		ϕᵣ₊₂ = ψ.nb[v_pos,h_pos,z_pos].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[v_pos,h_pos,z_pos].ϕᵣ₊₃
        h = fᵣ(ϕ,ψ.nb[v_pos,h_pos,z_pos],h_pos,ψ.consts)              # Higgs terms
#		@test h == -0.5*(ϕ.u⁺^2+ϕ.u⁺^2+ϕ.u⁻^2+ϕ.u⁻^2) + (ϕ.u⁺*ϕ.u⁻)^2*(2+ψ.consts.ν) + 0.5*(ϕ.u⁺^4 + ϕ.u⁻^4)
		energy += h
        energy += maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, g⁻²)
    end
    energy
end

# --------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
# This is general ordinary lattice regularization energy without CP1.
function ΔE(ϕ′::LatticeSite, ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, nnnb::NNNeighbors, h_pos::Int64, c::SystConstants)
    δE::Float64 = 0.0

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
        A⁰₋ = two_pi*c.f*(c.L-1)
    else
        A⁰₋ = two_pi*c.f*(h_pos-2)
    end
    
    # Normal kinetic terms
    δFₖ = 2*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁺-ϕ′.A[1]) - ϕᵣ₋₁.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-ϕ.A[1]) - ϕᵣ₋₁.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]))
            + (ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁺-ϕ′.A[2]-A⁰) - ϕᵣ₋₂.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰)) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-ϕ.A[2]-A⁰) - ϕᵣ₋₂.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-ϕᵣ₋₂.A[2]-A⁰))
     + c.κ₅*((ϕ′.u⁺^2 - ϕ′.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-ϕ′.A[3]) - ϕᵣ₋₃.u⁺*ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3])) 
             - (ϕ.u⁺^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3]) - ϕᵣ₋₃.u⁺*ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3])))
           + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁻-ϕ′.A[1]) - ϕᵣ₋₁.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1])) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-ϕ.A[1]) - ϕᵣ₋₁.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
            + (ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁻-ϕ′.A[2]-A⁰) - ϕᵣ₋₂.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A[2]-A⁰)) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-ϕ.A[2]-A⁰) - ϕᵣ₋₂.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-ϕᵣ₋₂.A[2]-A⁰))
     + c.κ₅*((ϕ′.u⁻^2 - ϕ′.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ′.θ⁻-ϕ′.A[3]) - ϕᵣ₋₃.u⁻*ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3])) 
             - (ϕ.u⁻^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A[3]) - ϕᵣ₋₃.u⁻*ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3]))))

    # Potential terms
    δFᵥ = (((ϕ′.u⁺*ϕ′.u⁻)^2*(2+c.ν*cos(2*(ϕ′.θ⁺-ϕ′.θ⁻))) - (1-1/c.β)*ϕ′.u⁺^2 + 0.5*ϕ′.u⁺^4 - (1-1/c.β)*ϕ′.u⁻^2 + 0.5*ϕ′.u⁻^4)
           - ((ϕ.u⁺*ϕ.u⁻)^2*(2+c.ν*cos(2*(ϕ.θ⁺-ϕ.θ⁻))) - (1-1/c.β)*ϕ.u⁺^2 + 0.5*ϕ.u⁺^4 - (1-1/c.β)*ϕ.u⁻^2 + 0.5*ϕ.u⁻^4))

    # Andreev-Bashkin terms
    δFₐₙ = c.γ^2*(c.ν+1)*(ϕᵣ₊₂.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-(ϕ′.A[2]+A⁰)) - ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-(ϕ.A[2]+A⁰)))
        - ϕᵣ₊₁.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-ϕ′.A[1]) - ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₂.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-(ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-(ϕ.A[2]+A⁰))) 
        - ϕᵣ₊₁.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])))
    
    # Mixed gradient terms
    δFₘ = c.γ^2*(c.ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*(sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ - (ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
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
    #δE += ((ϕ′.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ′.A[2])^2 - (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])^2)*c.g⁻²
    δFₐ = maxwell(ϕ′, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²) - maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²)
    # Then from position r-x
    δFₐ += maxwell(ϕᵣ₋₁, ϕ′, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²) - maxwell(ϕᵣ₋₁, ϕ, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²)
    #δE += ((ϕᵣ₋₁.A[1] + ϕ′.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2 - (ϕᵣ₋₁.A[1] + ϕ.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2)*c.g⁻²
    # Then from position r-y
    δFₐ += maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ′, ϕᵣ₋₂₊₃, c.g⁻²) - maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ, ϕᵣ₋₂₊₃, c.g⁻²)
    #δE += ((ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ′.A[1] - ϕᵣ₋₂.A[2])^2 - (ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ.A[1] - ϕᵣ₋₂.A[2])^2)*c.g⁻²
    # Then we also have to include the position r-z
    δFₐ += maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ′, c.g⁻²) - maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ, c.g⁻²)

	δE = δFₖ + δFᵥ + δFₐ# + δFₐₙ + δFₘ 
end

# --------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
function ΔECP1(ϕ′::LatticeSite, ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, nnnb::NNNeighbors, h_pos::Int64, c::SystConstants)
    δE::Float64 = 0.0

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
        A⁰₋ = two_pi*c.f*(c.L-1)
    else
        A⁰₋ = two_pi*c.f*(h_pos-2)
    end
    
    # Normal kinetic terms
    δE += -2*c.γ^2*(
		  ϕᵣ₊₁.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ′.θ⁺ - ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - ϕ.A[1]))
        + ϕᵣ₋₁.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₁.u⁻*(ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ′.θ⁻ - ϕ′.A[1]) - ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - ϕ.A[1]))
        + ϕᵣ₋₁.u⁻*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))

        + ϕᵣ₊₂.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ′.θ⁺ - (ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - (ϕ.A[2]+A⁰)))
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)))
        + ϕᵣ₊₂.u⁻*(ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ′.θ⁻ - (ϕ′.A[2]+A⁰)) - ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - (ϕ.A[2]+A⁰)))
		+ ϕᵣ₋₂.u⁻*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰))))

    # Normal z-direction kinetic term
    δE += 2*c.κ₅*(
                  ϕ′.u⁺^2 - ϕ.u⁺^2
                + ϕ′.u⁻^2 - ϕ.u⁻^2
                - ϕᵣ₊₃.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ′.θ⁺-ϕ′.A[3]) - ϕ.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-ϕ.A[3]))
                - ϕᵣ₊₃.u⁻*(ϕ′.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ′.θ⁻-ϕ′.A[3]) - ϕ.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-ϕ.A[3]))
                - ϕᵣ₋₃.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₃.θ⁺-ϕᵣ₋₃.A[3]))
                - ϕᵣ₋₃.u⁻*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₃.θ⁻-ϕᵣ₋₃.A[3])))
    
    # Potential energy terms
    δE += c.γ^4*((ϕ′.u⁺*ϕ′.u⁻)^2*(1+c.ν*cos(2*(ϕ′.θ⁺ - ϕ′.θ⁻))) - (ϕ.u⁺*ϕ.u⁻)^2*(1+c.ν*cos(2*(ϕ.θ⁺ - ϕ.θ⁻))))
    
    # Andreev-Bashkin terms
    δE += c.γ^2*(c.ν+1)*(ϕᵣ₊₂.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻-(ϕ′.A[2]+A⁰)) - ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-(ϕ.A[2]+A⁰)))
        - ϕᵣ₊₁.u⁺*(ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻-ϕ′.A[1]) - ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁻*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₂.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺-(ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-(ϕ.A[2]+A⁰))) 
        - ϕᵣ₊₁.u⁻*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺-ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-ϕ.A[1])) 
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁺*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1])))
    
    # Mixed gradient terms
    δE += c.γ^2*(c.ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*(sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ - (ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
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
    #δE += ((ϕ′.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ′.A[2])^2 - (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])^2)*c.g⁻²
    δE += maxwell(ϕ′, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²) - maxwell(ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₃, c.g⁻²)
    # Then from position r-x
    δE += maxwell(ϕᵣ₋₁, ϕ′, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²) - maxwell(ϕᵣ₋₁, ϕ, ϕᵣ₋₁₊₂, ϕᵣ₋₁₊₃, c.g⁻²)
    #δE += ((ϕᵣ₋₁.A[1] + ϕ′.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2 - (ϕᵣ₋₁.A[1] + ϕ.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2)*c.g⁻²
    # Then from position r-y
    δE += maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ′, ϕᵣ₋₂₊₃, c.g⁻²) - maxwell(ϕᵣ₋₂, ϕᵣ₊₁₋₂, ϕ, ϕᵣ₋₂₊₃, c.g⁻²)
    #δE += ((ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ′.A[1] - ϕᵣ₋₂.A[2])^2 - (ϕᵣ₋₂.A[1] + ϕᵣ₊₁₋₂.A[2] - ϕ.A[1] - ϕᵣ₋₂.A[2])^2)*c.g⁻²
    # Then we also have to include the position r-z
    δE += maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ′, c.g⁻²) - maxwell(ϕᵣ₋₃, ϕᵣ₊₁₋₃, ϕᵣ₊₂₋₃, ϕ, c.g⁻²)

	δE
end

# -------------------------------------------------------------------------------------------------
# Function for testing energy difference, made for testing E vs ΔE functions.
function testEnDiff(ψ₀::State, ψ::State, dE::Float64)
    E₀ = E(ψ₀)
    E′ = E(ψ)
    if E₀ > E′
        same = isapprox(E′-dE, E₀; atol=0, rtol=1e-13)
    else
        same = isapprox(E′, E₀+dE; atol=0, rtol=1e-13)
    end
    if !same
        println("|(E′-E₀)-dE|/max(E′-E₀, dE) = $(abs((E′-E₀)-dE)/max(E′-E₀, dE))")
        println("E₀\t\t\tΔE\t\t\tE′-E₀\t\t\t|ΔE-E₁-E₂|")
        println("$E₀\t$dE\t$(E′-E₀)\t$(abs(dE - (E′-E₀)))")
    end
    return same
end

