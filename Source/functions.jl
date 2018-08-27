using Distributions
using StatsBase

####################################################################################################################
#                            Misc. help functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Helping function for calculating position vector in [x,y] format from a lattice position
# assuming origo is in lower left corner of the square lattice. L is the size of the lattice along one dimension.
function getVectorPosition(L::Int64, pos::Array{Int64,1})
    v_pos = pos[1]
    h_pos = pos[2]
    return [h_pos-1, L-v_pos]
end

# -----------------------------------------------------------------------------------------------------------
function autocorrTime(O::Array{Float64,1}, c::Float64=5.0)
    N = size(O,1)
    
    # Estimate all correlation functions for different lags=[0:N]
    ρ = autocor(O, 0:(N-1))
    # Estimate correlation times. Each element in the array is an estimate for the correlation time τ(m) 
    # where summing over the i=1:M lagged autocorrelation functions using different numbers  M <= N
    τ = 2*cumsum(ρ)-1
    
    for m=1:N
        if m < c*τ[m]
            # We return the estimate for τ that is such that M < c*τ(m) 
            return τ[m]
        end
    end
    return τ[N]
end

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


####################################################################################################################
#                            Monte-Carlo functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Given a lattice site ϕ, propose a new lattice site with values in intervals around the existing ones.
function proposeLocalUpdate(ϕ::LatticeSite, sim::Controls)
    
    u⁺ = mod(ϕ.u⁺ + rand(Uniform(-sim.umax,sim.umax)),1) # This does not allow u⁺ = 1, is this a problem?
    # Construct new configuration at lattice site.
    return LatticeSite([ϕ.A[1]+rand(Uniform(-sim.Amax,sim.Amax)), ϕ.A[2]+rand(Uniform(-sim.Amax,sim.Amax))],
        mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
        u⁺, √(1-u⁺^2))
end

# -----------------------------------------------------------------------------------------------------------
# Performes a Metropolis Hasting update on a lattice site at position pos in state ψ given an inverse temperature
# β and where ϕᵣ... gives nearest and next nearest neighbor sites. Note that pos gives [y,x] of the position of
# the lattice site in normal array notation such that [1,1] is the upper left corner.
function metropolisHastingUpdate!(ψ::State, pos::Array{Int64,1}, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite,
        ϕᵣ₋₁::LatticeSite, ϕᵣ₋₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, ϕᵣ₋₂₊₁::LatticeSite, sim::Controls)
	# Save the lattice site at the targeted position in a temporary variable ϕ and use the lattice site
	# as a basis for proposing a new lattice site ϕ′. Then find the energy difference between having
	# ϕ′ or ϕ at position pos.
    const ϕ = ψ.lattice[pos...]
    const ϕ′ = proposeLocalUpdate(ϕ, sim)
    const δE = ΔE(ϕ′, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, pos[2], ψ.consts)
    
    # Create random number ran ∈ (0,1].
    const ran = rand()
    if ran==0
        ran=1
    end
    
    # Update state with probability min(1, e^{-β⋅δE})
    # and return the energy of final state regardless of whether it gets updated or not.
    if log(ran) <= -ψ.consts.β*δE
        ψ.lattice[pos...] = ϕ′
        return δE
    else
        return 0.0
    end
end

# -----------------------------------------------------------------------------------------------------------
# Takes a state ψ with an L×L lattice and tries to update each site on the lattice by running the
# metropolisHastingUpdate! function on it. Each part of the boundary is updated separately so that periodic
# boundary conditions are taken care of for values stored in each lattice site.
function mcSweep!(ψ::State, sim::Controls = Controls(π/3, 0.4, 3.0))
   
    # Find size of the lattice L
    const L::Int64 = ψ.consts.L
    
    
    # Updating upper right corner
    const ϕᵣ₊₁ = ψ.lattice[1,1]
    const ϕᵣ₊₂ = ψ.lattice[L,L]
    const ϕᵣ₋₁ = ψ.lattice[1,L-1]
    const ϕᵣ₋₂ = ψ.lattice[2,L]
    const ϕᵣ₋₁₊₂ = ψ.lattice[L,L-1]
    const ϕᵣ₋₂₊₁ = ψ.lattice[2,1]
    metropolisHastingUpdate!(ψ, [1,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating boundary paralell to y-axis
    # except for the upper and lower right corner.
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,1]
        ϕᵣ₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₁ = ψ.lattice[y,L-1]
        ϕᵣ₋₂ = ψ.lattice[y+1,L]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,1]
        metropolisHastingUpdate!(ψ, [y,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating lower right corner
    ϕᵣ₊₁ = ψ.lattice[L,1]
    ϕᵣ₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₁ = ψ.lattice[L,L-1]
    ϕᵣ₋₂ = ψ.lattice[1,L]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L-1]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,1]
    metropolisHastingUpdate!(ψ, [L,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating lower x-axis boundary except for lower left and right corner
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[L,x+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,x]
        ϕᵣ₋₁ = ψ.lattice[L,x-1]
        ϕᵣ₋₂ = ψ.lattice[1,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L-1,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[1,x+1]
        metropolisHastingUpdate!(ψ, [L,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating lower left corner
    ϕᵣ₊₁ = ψ.lattice[L,2]
    ϕᵣ₊₂ = ψ.lattice[L-1,1]
    ϕᵣ₋₁ = ψ.lattice[L,L]
    ϕᵣ₋₂ = ψ.lattice[1,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,2]
    metropolisHastingUpdate!(ψ, [L,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating left y-axis boundary except corners
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,2]
        ϕᵣ₊₂ = ψ.lattice[y-1,1]
        ϕᵣ₋₁ = ψ.lattice[y,L]
        ϕᵣ₋₂ = ψ.lattice[y+1,1]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,2]
        metropolisHastingUpdate!(ψ, [y,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating upper left corner
    ϕᵣ₊₁ = ψ.lattice[1,2]
    ϕᵣ₊₂ = ψ.lattice[L,1]
    ϕᵣ₋₁ = ψ.lattice[1,L]
    ϕᵣ₋₂ = ψ.lattice[2,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[2,2]
    metropolisHastingUpdate!(ψ, [1,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating upper x-axis boundary
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[1,x+1]
        ϕᵣ₊₂ = ψ.lattice[L,x]
        ϕᵣ₋₁ = ψ.lattice[1,x-1]
        ϕᵣ₋₂ = ψ.lattice[2,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[2,x+1]
        metropolisHastingUpdate!(ψ, [1,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # This concludes the updates of the boundary
    
    # Update the bulk
    for x=2:(L-1), y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,x+1]
        ϕᵣ₊₂ = ψ.lattice[y-1,x]
        ϕᵣ₋₁ = ψ.lattice[y,x-1]
        ϕᵣ₋₂ = ψ.lattice[y+1,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,x+1]
        metropolisHastingUpdate!(ψ, [y,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
end

# -----------------------------------------------------------------------------------------------------------
# Same as above but returns the fraction of accepted proposals over number of proposals and calculates
# nearest neighbors in a dynamic but more costly way.
function mcSweepFrac!(ψ::State, sim::Controls = Controls(π/3, 0.4, 1.0))
    count = 0
    L = ψ.consts.L
    for h_pos = 1:L
        for v_pos = 1:L
            ϕ = ψ.lattice[v_pos, h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,mod(h_pos,L)+1]
            ϕᵣ₊₂ = ψ.lattice[mod(v_pos-2,L)+1,h_pos]
            ϕᵣ₋₁ = ψ.lattice[v_pos,mod(h_pos-2,L)+1]
            ϕᵣ₋₂ = ψ.lattice[mod(v_pos,L)+1,h_pos]
            ϕᵣ₋₁₊₂ = ψ.lattice[mod(v_pos-2,L)+1, mod(h_pos-2,L)+1]
            ϕᵣ₋₂₊₁ = ψ.lattice[mod(v_pos,L)+1, mod(h_pos,L)+1]

            if metropolisHastingUpdate!(ψ,[v_pos,h_pos],ϕᵣ₊₁,ϕᵣ₊₂,ϕᵣ₋₁,ϕᵣ₋₂,ϕᵣ₋₁₊₂,ϕᵣ₋₂₊₁, sim) != 0.0
                count += 1
            end
        end
    end
    return count/L^2
end


# -----------------------------------------------------------------------------------------------------------
# Same as mcSweep but adds up the energy-differences and returns it.
function mcSweepEn!(ψ::State, sim::Controls = Controls(π/3, 0.4, 3.0))
    
    δE = 0.0
   
    # Find size of the lattice L
    const L::Int64 = ψ.consts.L
    
    
    # Updating upper right corner
    const ϕᵣ₊₁ = ψ.lattice[1,1]
    const ϕᵣ₊₂ = ψ.lattice[L,L]
    const ϕᵣ₋₁ = ψ.lattice[1,L-1]
    const ϕᵣ₋₂ = ψ.lattice[2,L]
    const ϕᵣ₋₁₊₂ = ψ.lattice[L,L-1]
    const ϕᵣ₋₂₊₁ = ψ.lattice[2,1]
    δE += metropolisHastingUpdate!(ψ, [1,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating boundary paralell to y-axis
    # except for the upper and lower right corner.
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,1]
        ϕᵣ₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₁ = ψ.lattice[y,L-1]
        ϕᵣ₋₂ = ψ.lattice[y+1,L]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,1]
        δE += metropolisHastingUpdate!(ψ, [y,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating lower right corner
    ϕᵣ₊₁ = ψ.lattice[L,1]
    ϕᵣ₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₁ = ψ.lattice[L,L-1]
    ϕᵣ₋₂ = ψ.lattice[1,L]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L-1]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,1]
    δE += metropolisHastingUpdate!(ψ, [L,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating lower x-axis boundary except for lower left and right corner
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[L,x+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,x]
        ϕᵣ₋₁ = ψ.lattice[L,x-1]
        ϕᵣ₋₂ = ψ.lattice[1,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L-1,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[1,x+1]
        δE += metropolisHastingUpdate!(ψ, [L,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating lower left corner
    ϕᵣ₊₁ = ψ.lattice[L,2]
    ϕᵣ₊₂ = ψ.lattice[L-1,1]
    ϕᵣ₋₁ = ψ.lattice[L,L]
    ϕᵣ₋₂ = ψ.lattice[1,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,2]
    δE += metropolisHastingUpdate!(ψ, [L,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating left y-axis boundary except corners
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,2]
        ϕᵣ₊₂ = ψ.lattice[y-1,1]
        ϕᵣ₋₁ = ψ.lattice[y,L]
        ϕᵣ₋₂ = ψ.lattice[y+1,1]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,2]
        δE += metropolisHastingUpdate!(ψ, [y,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # Updating upper left corner
    ϕᵣ₊₁ = ψ.lattice[1,2]
    ϕᵣ₊₂ = ψ.lattice[L,1]
    ϕᵣ₋₁ = ψ.lattice[1,L]
    ϕᵣ₋₂ = ψ.lattice[2,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[2,2]
    δE += metropolisHastingUpdate!(ψ, [1,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    
    # Updating upper x-axis boundary
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[1,x+1]
        ϕᵣ₊₂ = ψ.lattice[L,x]
        ϕᵣ₋₁ = ψ.lattice[1,x-1]
        ϕᵣ₋₂ = ψ.lattice[2,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[2,x+1]
        δE += metropolisHastingUpdate!(ψ, [1,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    
    # This concludes the updates of the boundary
    
    # Update the bulk
    for x=2:(L-1), y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,x+1]
        ϕᵣ₊₂ = ψ.lattice[y-1,x]
        ϕᵣ₋₁ = ψ.lattice[y,x-1]
        ϕᵣ₋₂ = ψ.lattice[y+1,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,x+1]
        δE += metropolisHastingUpdate!(ψ, [y,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, sim)
    end
    return δE
end


# -------------------------------------------------------------------------------------------------
# Takes a state ψ and plots the fraction of accepted proposals on number of proposals for M
# Monte Carlo sweeps of the lattice. Then return the average and standard deviation of this
# fraction as well as the time-series itself.
function mcProposalFraction(ψ::State, sim::Controls=Controls(π/2, 0.4, 3.0), M::Int64=500)
	ψ_copy = copy(ψ)
    L = ψ.consts.L
    fracs = zeros(M)
    
    # Go through the entire lattice M times and gain the statistic of whether it gets updated or not
    for i = 1:M
        fracs[i] = mcSweepFrac!(ψ_copy, sim)
    end
    res = mean(fracs)
    stdev = std(fracs)
    return (res, stdev, fracs)
end

# Same as mcProposalFraction but does update the state.
# The number of times the state has been updated can be surmised from the input parameter M.
# Does not return the series of acceptance-rates unlike mcProposalFraction
function mcProposalFraction!(ψ::State, sim::Controls=Controls(π/2, 0.4, 3.0), M::Int64=500)
    L = ψ.consts.L
    fracs = zeros(M)
    
    # Go through the entire lattice M times and gain the statistic of whether it gets updated or not
    for i = 1:M
        fracs[i] = mcSweepFrac!(ψ, sim)
    end
    av = mean(fracs)
    stdev = std(fracs)
    return (av, stdev)
end


# -----------------------------------------------------------------------------------------------------------
# Adjusts the simulation controls such that ψ has an AR larger than LOWER. Tries to find the simulation controls
# that are simultaneously the ones with the highest values.
function adjustSimConstants!(sim::Controls, ψ::State, M::Int64 = 40)
    #println("Adjusting simulation constants $(sim.θmax), $(sim.umax), $(sim.Amax)")
    CUTOFF_MAX = 42       # How many times the while loop should run. 
    LOWER = 0.3           # Minimum acceptance rate (AR)
    NEEDED_PROPOSALS = 5  # Number of sim with AR >= LOWER
    DIVI = 1.5            # The number the current value is divided by to get lower limit of search interval.
    TRIED_VALUES = 4      # Number of values to try in new interval
    
    proposedConstants = [copy(sim) for i=1:NEEDED_PROPOSALS]
    proposedAR = zeros(NEEDED_PROPOSALS)
    proposals = 0
    n = 0
    adjustment_mcs = 0
    
    # First we get an estimate of the accept probability
    (av, std) = mcProposalFraction!(ψ, sim, M)
    adjustment_mcs = M
    if av >= LOWER
        return (av, adjustment_mcs)
    end
    
    s₀ = copy(sim)
    while (proposals < NEEDED_PROPOSALS)
        # The starting state of this loop will be that we have no proposals for sims that has acceptance rate higher
        # than LOWER (including the initial s₀)
        
        # First we look at Amax
        interval_end = s₀.Amax/DIVI
        tries_A = [Controls(s₀.θmax, s₀.umax, 
                s₀.Amax-x*(s₀.Amax-interval_end)/TRIED_VALUES) for x = 1:TRIED_VALUES]
        f_A = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_A[i], std) = mcProposalFraction!(ψ, tries_A[i], M)
            adjustment_mcs += M
            # If we find a value with probability >= LOWER, then this is the largest such value we have found
            # and should be included in proposed constants
            if f_A[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_A[i]
                proposedAR[proposals] = f_A[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break

        # Then we try to vary umax
        interval_end = s₀.umax/DIVI
        tries_u = [Controls(s₀.θmax, s₀.umax - x*(s₀.umax-interval_end)/TRIED_VALUES,
                s₀.Amax) for x = 1:TRIED_VALUES]
        f_u = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_u[i], std) = mcProposalFraction!(ψ, tries_u[i], M)
            adjustment_mcs += M
            if f_u[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_u[i]
                proposedAR[proposals] = f_u[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break

        # Then we try to vary θmax
        interval_end = s₀.θmax/DIVI
        tries_θ = [Controls(s₀.θmax - x*(s₀.θmax - interval_end)/TRIED_VALUES, s₀.umax,
                s₀.Amax) for x = 1:TRIED_VALUES]
        f_θ = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_θ[i], std) = mcProposalFraction!(ψ, tries_θ[i], M)
            adjustment_mcs += M
            if f_θ[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_θ[i]
                proposedAR[proposals] = f_θ[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break
        
        # At this point in the loop we have tried to get some proposals but failed to get enough of them.
        # We need to start the loop again with a new starting state that is such that the accept ratio is
        # as high is possible so that we can get more proposals.
        accept_ratios = vcat(f_A, f_u, f_θ)
        tries = vcat(tries_A, tries_u, tries_θ)
        i_max = indmax(accept_ratios)   # Finding index of sim that gave highest accept ratio.
        s₀ = tries[i_max]               # Setting this sim to the initial one.
        
        
        n += 1
        if n >= CUTOFF_MAX
            println("WARNING: Could not find simulation constant such that update probability 
                was higher than $(LOWER)")
            sim = s₀
            return (accept_ratios[i_max], adjustment_mcs)
        end
    end
    # The end situation of the loop is that we have a number of proposals >= NEEDED_PROPOSALS.
    
    # Finding the distance from zero of the different proposals.
    norms = zeros(proposals)
    for i = 1:proposals
        norms[i] = proposedConstants[i].θmax^2 + proposedConstants[i].umax^2 + proposedConstants[i].Amax^2
    end
    i_max = indmax(norms)
    setValues!(sim, proposedConstants[i_max]) # Finally we update the simulation constants to the sim that has 
    # highest norm, and an accept ratio above LOWER.
    return (proposedAR[i_max], adjustment_mcs)
    # Return the acceptance ratio of the new sim and the number of Monte-Carlo Sweeps done during this adjustment.
end


# -----------------------------------------------------------------------------------------------------------
# Version of find Equilibrium that doesn't waste MCS when adjusting sim Constants and also calculates energy based
# on return value of mcSweepEn! as well as gets rid of the use of dynamic arrays.
# -----------------------------------------------------------------------------------------------------------
# Given values for the physical constants of the system as well as the system size, we find the number of MC-sweeps it
# takes until the internal energy of the system reaches a more or less constant value.
function findEquilibrium(c::SystConstants, sim₁::Controls=Controls(π/3, 0.4, 3.0), 
        T::Int64=1000, ex::Float64=1.5, di::Int64=8)
    CUTOFF_MAX::Int64 = 40000
    ADJUST_INTERVAL::Int64 = 2000
    STD_NUMBER::Int64 = 1
    println("Finding Equilibrium of\n$(c)\n$(sim₁)")
    
    ψ₂ = State(2, c)
    sim₂ = copy(sim₁)
    ψ₁ = State(1, c)
    dE = zeros(CUTOFF_MAX)
    E₁ = zeros(CUTOFF_MAX)
    E₂ = zeros(CUTOFF_MAX)
    adjustment_mcs = 0
    
    # Check that the un-correllated state has higher energy than the correlated
    if E(ψ₂) <= E(ψ₁)
        error("Correlated state has higher energy than un-correlated")
    end
    
    tₛ = 0    # The wanted t₀ does not exist at or before this position.
    t₀ = T
    
    E₁[1] = E(ψ₁)
    E₂[1] = E(ψ₂)
    dE[1] = E₂[1] - E₁[1]
    for i = 2:T
        E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
        E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
        dE[i] = E₂[i] - E₁[i]
    end
    
    # Adjust simulation constants as needed
    (ar, mcs1) = adjustSimConstants!(sim₁, ψ₁)
    (ar, mcs2) = adjustSimConstants!(sim₂, ψ₂)
    adjustment_mcs += max(mcs1, mcs2)   # Adds the max number of monte-carlo sweeps done for the two
    # adjustments to the number of adjustments used for finding the equilibrium time.
    
    E₁[T] = E(ψ₁)
    E₂[T] = E(ψ₂)
    dE[T] = E₂[T] - E₁[T]
    
    while tₛ < CUTOFF_MAX
        # Find the first occurence of dE <= 0 if it exists
        println("Searching for ΔE <= 0..")
        t₀ = T
        for i = (tₛ+1):T
            if dE[i] <= 0
                t₀ = i
                break
            end
        end
        
        while T <= t₀ && T < CUTOFF_MAX
            # If we couldn't find a t₀ in dE we have to try and increase simulation time
            tₛ = T
            T = min(Int(ceil(T*ex)), CUTOFF_MAX)
            for i = (tₛ+1):T
                print("$(round(i/CUTOFF_MAX*100,1))% of max\r") # Debug
                E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
                E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
                dE[i] = E₂[i] - E₁[i]
                
                # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
                if i % ADJUST_INTERVAL == 0
                    ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                    ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                    adjustment_mcs += max(mcs1, mcs2)
                    
                    # Then we also go over an extra time to get energy correct
                    E₁[i] = E(ψ₁)
                    E₂[i] = E(ψ₂)
                    dE[i] = E₂[i] - E₁[i]
                end
            end
            
            # Then we again see if we can find the first occurrence of dE <= 0 after tₛ
            t₀ = T
            for i = tₛ:T
                if dE[i] <= 0
                    t₀ = i
                    break
                end
            end
            
            if t₀ == T == CUTOFF_MAX # We have not found any dE < 0 and we have reached the max number of sweeps
                println("Failed to find a point where ΔE <= 0")
                return (-1, E₁, E₂, dE, ψ₁, ψ₂, sim₁, sim₂)
            end
            
            # When the loop ends we should have the situation that T > t₀ where t₀ is the first occurrence
            # in dE where dE[t₀] <= 0
        end
        println("ΔE <= 0 found at t₀ = $(t₀)!\nChecking if average is close to 0..")
        
        # Now we make sure that T is large enough such that [1,T] includes an interval [t₀, t₀+t₀/div]
        # so that an average can be performed
        t_end = min(t₀ + Int(ceil(t₀/di)), CUTOFF_MAX)
        while T < t_end
            T += 1
            E₁[T] = E₁[T-1] + mcSweepEn!(ψ₁, sim₁)
            E₂[T] = E₂[T-1] + mcSweepEn!(ψ₂, sim₂)
            dE[T] = E₂[T] - E₁[T]

            # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
            if T % ADJUST_INTERVAL == 0
                ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                adjustment_mcs += max(mcs1, mcs2)
                
                # Then we also go over an extra time to get energy correct
                E₁[i] = E(ψ₁)
                E₂[i] = E(ψ₂)
                dE[i] = E₂[i] - E₁[i]
            end
        end
        
        # Now we calculate the average and standard deviation of dE over [t₀, T] and check if the
        # average is within STD_NUMBER of standard deviations of 0 at which point we declare the equilibrium
        # time to have been found.
        int = dE[t₀:T]
        av = mean(int)
        st = std(int)
        if abs(av) <= STD_NUMBER*st
            println("Equilibrium found at time $(T+adjustment_mcs)
over the interval [$(t₀), $(T)]
s.t. <ΔE> = $(round(av,2)) ± $(round(st/sqrt(size(int,1)), 1))
std(ΔE) = $(round(st, 1))")
            return (T+adjustment_mcs, E₁[1:T], E₂[1:T], dE[1:T], ψ₁, ψ₂, sim₁, sim₂)
        end
        
        println("Average was not close to 0. Increasing interval.")
        
        # If we didn't find an interval that had an average close to 0 we assume this interval is ahead of us
        # and start again with an increased T, setting the starting point tₛ to the end of the interval.
        tₛ = T
        T = min(Int(ceil(T*ex)), CUTOFF_MAX)
        
        # Simulating new MCS
        for i = (tₛ+1):T
            E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
            E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
            dE[i] = E₂[i] - E₁[i]

            # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
            if i % ADJUST_INTERVAL == 0
                ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                adjustment_mcs += max(mcs1, mcs2)
                
                # Then we also go over an extra time to get energy correct
                E₁[i] = E(ψ₁)
                E₂[i] = E(ψ₂)
                dE[i] = E₂[i] - E₁[i]
            end
        end
    end
    return (-1, E₁, E₂, dE, ψ₁, ψ₂, sim₁, sim₂)
end
# IDEA: We could have sim be controlled as a Metropolis Update where the internal energy equivalent could be the number of
# accepts pr proposal over mcSweep and how close that is to 1/2.




####################################################################################################################
#                            Planar Structure function
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Returns the un-normalized local vorticity by preforming a plaquette sum using the gauge-invariant
# difference of the θ field.
function n⁺(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    return (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos) 
        - mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺, two_pi) + ϕᵣ₊₂.A[1] 
        - mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺, two_pi) + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
end
function n⁻(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    return (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos) 
        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻, two_pi) + ϕᵣ₊₂.A[1] 
        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻, two_pi) + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
end

# -----------------------------------------------------------------------------------------------------------
function structureFunctionPluss{T<:Real}(k::Array{T,1}, ψ::State)
    sum = Complex(0)
    L = ψ.consts.L
    
    # Sum over the corners
    # Upper left corner
     r = [0, L-1] # For r we assume origo is in position [L,1] of the lattice. 
                  # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
                  # Additionally we define it such that we get the usual r = [x,y] order of dimensions.
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower left corner
     r = [0,0]
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower right corner
     r = [L-1,0]
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Upper right corner
     r = [L-1,L-1]
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        r = [i-1, L-1]
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Left border (counting in -y direction)
        r = [0, L-i]
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        
        # Lower border (counting in +x direction)
        r = [i-1, 0]
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Right border (counting in -y direction)
        r = [L-1, L-i]
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    
    return abs2(sum)
end
function structureFunctionMinus{T<:Real}(k::Array{T,1}, ψ::State)
    sum = Complex(0)
    L = ψ.consts.L
    
    # Sum over the corners
    # Upper left corner
     r = [0, L-1] # For r we assume origo is in position [L,1] of the lattice. 
                  # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
                  # Additionally we define it such that we get the usual r = [x,y] order of dimensions.
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower left corner
     r = [0,0]
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower right corner
     r = [L-1,0]
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Upper right corner
     r = [L-1,L-1]
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        r = [i-1, L-1]
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Left border (counting in -y direction)
        r = [0, L-i]
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        
        # Lower border (counting in +x direction)
        r = [i-1, 0]
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Right border (counting in -y direction)
        r = [L-1, L-i]
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            sum += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    
    return abs2(sum)
end
function structureFunction{T<:Real}(k::Array{T,1}, ψ::State)
    sum⁺ = Complex(0)
    sum⁻ = Complex(0)
    L = ψ.consts.L
    
    # Sum over the corners
    # Upper left corner
     r = [0, L-1] # For r we assume origo is in position [L,1] of the lattice. 
                  # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
                  # Additionally we define it such that we get the usual r = [x,y] order of dimensions.
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower left corner
     r = [0,0]
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower right corner
     r = [L-1,0]
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Upper right corner
     r = [L-1,L-1]
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        r = [i-1, L-1]
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Left border (counting in -y direction)
        r = [0, L-i]
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        
        # Lower border (counting in +x direction)
        r = [i-1, 0]
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Right border (counting in -y direction)
        r = [L-1, L-i]
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
            sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    
    return (abs2(sum⁺),abs2(sum⁻))
end

# -----------------------------------------------------------------------------------------------------------
# This version should do the same as above, but determines the neighbors dynamically at a cost to performance.
function structureFunctionPlussDyn{T<:Real}(k::Array{T,1}, ψ::State)
    L = ψ.consts.L
    sum = Complex(0.0)
    
    # Sum over entire lattice and determine nearest neighbors dynamically.
    for h_pos = 1:L
        for v_pos = 1:L
            
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos, mod(h_pos, L)+1]
            ϕᵣ₊₂ = ψ.lattice[mod(v_pos-2,L)+1, h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[mod(v_pos-2,L)+1, mod(h_pos, L)+1]
            sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    return abs2(sum)
end


####################################################################################################################
#                            Thermal averages
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
function structureFunctionPlussAvg{T<:Real}(k::Array{T,1}, syst::SystConstants, M::Int64, Δt::Int64)
    S = zeros(M)
    secondMoment = 0
    s_norm = 1/(sim.f*sim.L^2*two_pi)^2
    # Finding a state that has reached thermal equilibrium
    (t₀, dE, ψ, ψ₂) = findEquilibrium(syst)
    
    # Loop over M measurements
    for m = 1:M
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ)
        end
        
        # Make a normalized measurement
        secondMoment = (S[m] = s_norm*structureFunctionPluss(k, ψ))^2
    end
    
    # Return average and error estimates
    av = mean(S)
    secondMoment = secondMoment/M
    τₒ = autocorrTime(S, 5.0)
    err = (1+2*τₒ)*(secondMoment - av^2)/(M-1)
    
    return (av, err, S)
end
# Take in a matrix of k-values and calculate both the vorticity of θ⁺ and θ⁻.
function structureFunctionAvg{T<:Real}(ks::Array{Array{T, 1}, 2}, syst::SystConstants, M::Int64, Δt::Int64)
    Lky = size(ks, 1)
    Lkx = size(ks, 2)
    S⁺ = [zeros(M) for y=1:Lky, x=1:Lkx]    # Matrix containing the series of measurements for each k
    S⁻ = [zeros(M) for y=1:Lky, x=1:Lkx]
    Sm⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    Sm⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    s_norm_inv = 1/(syst.f*syst.L^2*two_pi)^2
    
    # Finding a state that has eached thermal equilibrium
    println("Finding equilibrium")
    (t₀, dE, ψ, ψ₂) = findEquilibrium(syst)
    
    prinln("Making measurements")
    # Loop over M measurements
    for m = 1:M
        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ)
        end
        
        # Make measurements 
        for y = 1:Lky, x = 1:Lkx
            (S⁺[y,x][m], S⁻[y,x][m]) = structureFunction(ks[y,x], ψ)
            Sm⁺[y,x] += S⁺[y,x][m]^2
            Sm⁻[y,x] += S⁻[y,x][m]^2
        end
    end
    
    # Return average and error estimates
    avS⁺ = [mean(S⁺[y,x]) for y=1:Lky, x=1:Lkx]
    avS⁻ = [mean(S⁻[y,x]) for y=1:Lky, x=1:Lkx]
    τ⁺ = [autocorrTime(S⁺[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    τ⁻ = [autocorrTime(S⁻[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    errS⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    errS⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    
    for y=1:Lky, x=1:Lkx
        Sm⁺[y,x] /= M
        Sm⁻[y,x] /= M
        errS⁺[y,x] = (1+2*τ⁺[y,x])*(Sm⁺[y,x] - avS⁺[y,x]^2)/(M-1)
        errS⁻[y,x] = (1+2*τ⁻[y,x])*(Sm⁻[y,x] - avS⁻[y,x]^2)/(M-1)
    end
    
    return (avS⁺, errS⁺, avS⁻, errS⁻)
end


####################################################################################################################
#                            Diagnostic functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Testing how well our algorithms work for this choice of parameters.
function testSystem(syst::SystConstants, sim::Controls, M::Int64=2000)
    # Let's first look at the time it takes to equilibrate the state
    ψ = State(2, syst)
    ψ_old = copy(ψ)
    
    println("What happens after $(M) number of MCS")
    E_test = zeros(M)

    # Do mcSweep! M times
    for i = 1:M
        mcSweep!(ψ, sim)
        E_test[i] = E(ψ)
    end
    println("Checking if mcSweeped state has lower energy")
    println(@test E(ψ) < E(ψ_old))
    plt = plot(1:M, E_test, title="Development of internal free energy with Monte-Carlo sweeps", 
        xlabel="MCS", ylabel="E()")
    display(plt)
    
    # Then we find the proposal fraction for the update
    (res, stdev, fracs) = mcProposalFraction(ψ_old, M)
    plt = plot(1:M, fracs, title="Probability of update over time", ylabel="#accepted proposals / L^2", xlabel="MCS")
    display(plt)
    println("metropolisHastingUpdate! proposal was accepted $(round(res*100,1))±$(round(stdev*100,1))% of the times")
    
    # At infinite energy, metropolis Hasting will always accept the new state, therefore we expect
    # that after a couple of mcSweeps, the state will be completely changed
    println("Testing if mcSweep! gives completely different state when temperature is infinite, where completely different
    means that all values on all lattice sites are different. Thus we have proved that mcSweep! visits all lattice sites.")
    ψ = State(2, SystConstants(syst.L, syst.γ, syst.g⁻², syst.ν, syst.f, 0))
    ψ_old = copy(ψ)
    for i = 1:2
        mcSweep!(ψ)
    end
    println(@test latticeCompletelyDifferent(ψ,ψ_old))

    println("\nEquilibrium calculation\n----------------------------------------------------------------")
    ψ₁ = State(1, syst)
    ψ₂ = State(2, syst)
    println("Checking that a random state has lower energy than a completely correlated state")
    println(@test E(ψ₁) < E(ψ₂))
    flush(STDOUT)
    (t₀, dE, ψ₁, ψ₂) = findEquilibrium(syst, M)
    @show t₀
    T = size(dE,1)
    plt = plot(1:T,dE, title="Difference in energy for a correlated state - random state", xlabel="MCS", ylabel="E1-E2")
    display(plt)
end
