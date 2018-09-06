####################################################################################################################
#                            Energy functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms.
function fᵣ(ϕ::LatticeSite,ϕᵣ₊₁::LatticeSite,ϕᵣ₊₂::LatticeSite,A₀::Float64,c::SystConstants)
    energy = 0
    A₂ = ϕ.A[2]+A₀
    # Kinetic energy Fₖ
    energy += -2*c.γ^2*(ϕᵣ₊₁.u⁺ *ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺ - ϕ.A[1]) 
        + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺ - A₂) 
        + ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻ - ϕ.A[1]) 
        + ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻ - A₂) )
    # Potential energy Fᵥ
    energy += c.γ^4*ϕ.u⁺^2*ϕ.u⁻^2*(1+c.ν*cos(2*(ϕ.θ⁺-ϕ.θ⁻)))
    # Andreev Bashkin terms
    energy += c.γ^2*(c.ν+1)*(ϕ.u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ - A₂) 
        + ϕ.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ - A₂) 
        - ϕ.u⁻*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - ϕ.A[1]) 
        - ϕ.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.A[1]))
    # Mixed gradient terms
    energy += c.γ^2*(c.ν-1)*(ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*sin(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ - (ϕ.A[1]-A₂)) 
        - ϕᵣ₊₂.u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ - (ϕ.A[1] - A₂)) 
        + 2*ϕ.u⁺*ϕ.u⁻*sin(ϕ.θ⁻-ϕ.θ⁺) 
        + ϕ.u⁻*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ - A₂)    # Here there used to be a sign error
        - ϕ.u⁺*ϕᵣ₊₂.u⁻*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ - A₂)    # Here too
        + ϕ.u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ - ϕ.A[1]) 
        - ϕ.u⁺*ϕᵣ₊₁.u⁻*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.A[1]))
    energy
end

# -----------------------------------------------------------------------------------------------------------
# Loops over all positions of the lattice of a state and calculates the total energy from the
# Higgs-field terms using the function fᵣ() + the energy from the gauge field.
function E(ψ::State)
    const γ = ψ.consts.γ
    const g⁻² = ψ.consts.g⁻²
    energy = 0.0
    const ν = ψ.consts.ν
    const f = ψ.consts.f
    const L = ψ.consts.L
    
    # Contribution from upper right corner
    A⁰ = (L-1)*two_pi*f
    ϕ = ψ.lattice[1,L]    # Lattice site at upper right corner
    ϕᵣ₊₁ = ψ.lattice[1,1]    # Nearest neighbor at r+x is upper left corner
    ϕᵣ₊₂ = ψ.lattice[L,L]  # Nearest neighbor at r+y is lower right corner
    energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,ψ.consts)              # Higgs terms
    energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2*g⁻² # Maxwell term
    
    # Contribution from right boundary paralell to y-axis
    # except for the upper right corneϕ.
    for y=2:L
        ϕ = ψ.lattice[y,L]
        ϕᵣ₊₁ = ψ.lattice[y,1]
        ϕᵣ₊₂ = ψ.lattice[y-1,L]
        energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,ψ.consts)              # Higgs terms
        energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2*g⁻²
    end
    
    # Contribution from the bulk of lattice sites and upper boundary
    for x=1:(L-1)
        A⁰ = (x-1)*two_pi*f        # Constant vector potential.
        # Constribution from upper boundary except upper right corner
        ϕ = ψ.lattice[1,x]
        ϕᵣ₊₁ = ψ.lattice[1,x+1]
        ϕᵣ₊₂ = ψ.lattice[L,x]
        energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,ψ.consts)              # Higgs terms
        energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2*g⁻²
        
        # Contribution from the rest of the bulk.
        for y=2:L
            ϕ = ψ.lattice[y,x]          # Lattice site at position r
            ϕᵣ₊₁ = ψ.lattice[y,x+1]       # Nearest neighbor at r+x
            ϕᵣ₊₂ = ψ.lattice[y-1,x]       # Nearest neighbor at r+y
            
            energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,ψ.consts)              # Higgs terms
            energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2*g⁻²
        end
    end
    
    energy
end

# -----------------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
function ΔE(ϕ′::LatticeSite, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, 
        ϕᵣ₋₁::LatticeSite, ϕᵣ₋₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, ϕᵣ₋₂₊₁::LatticeSite, x::Int64, c::SystConstants)
    δE::Float64 = 0.0
    
    # Calculate constant link variables
    const A⁰ = two_pi*c.f*(x-1)
    #const A⁰₊ = 2π*c.f*x
    const A⁰₋ = two_pi*c.f*(x-2)
    
    # Normal kinetic terms
    δE += -2*c.γ^2*(ϕᵣ₊₁.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ′.θ⁺ - ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - ϕ.A[1]))
        + ϕᵣ₋₁.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₁.u⁻*(ϕ′.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ′.θ⁻ - ϕ′.A[1]) - ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - ϕ.A[1]))
        + ϕᵣ₋₁.u⁻*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻-ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₂.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ′.θ⁺ - (ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - (ϕ.A[2]+A⁰)))
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺-(ϕᵣ₋₂.A[2]+A⁰)))
        + ϕᵣ₊₂.u⁻*(ϕ′.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ′.θ⁻ - (ϕ′.A[2]+A⁰)) - ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - (ϕ.A[2]+A⁰)))
        + ϕᵣ₋₂.u⁻*(ϕ′.u⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻-(ϕᵣ₋₂.A[2]+A⁰))))
    
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
    # First contribution from current position
    δE += ((ϕ′.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ′.A[2])^2 - (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])^2)*c.g⁻²
    # Then from position r-x
    δE += ((ϕᵣ₋₁.A[1] + ϕ′.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2 - (ϕᵣ₋₁.A[1] + ϕ.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2)*c.g⁻²
    # Then from position r-y
    δE += ((ϕᵣ₋₂.A[1] + ϕᵣ₋₂₊₁.A[2] - ϕ′.A[1] - ϕᵣ₋₂.A[2])^2 - (ϕᵣ₋₂.A[1] + ϕᵣ₋₂₊₁.A[2] - ϕ.A[1] - ϕᵣ₋₂.A[2])^2)*c.g⁻²
end

