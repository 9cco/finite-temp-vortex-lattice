using Distributions

####################################################################################################################
#                            Energy functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Calculate energy contribution from a single term in the energy sum of the Higgs terms.
function fᵣ(ϕ::LatticeSite,ϕᵣ₊₁::LatticeSite,ϕᵣ₊₂::LatticeSite,A₀::Float64,γ::Float64,g::Float64,ν::Float64,f::Float64)
    energy = 0
    uᵣ₊₁⁻ = √(1-ϕᵣ₊₁.u⁺^2)
    u⁻ = √(1-ϕ.u⁺^2)
    uᵣ₊₂⁻ = √(1-ϕᵣ₊₂.u⁺^2)
    A₂ = ϕ.A[2]+A₀
    # Kinetic energy Fₖ
    energy += -2*γ^2*(ϕᵣ₊₁.u⁺ *ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺ + g*ϕ.A[1]) 
        + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺ + g*A₂) 
        + uᵣ₊₁⁻*u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻ + g*ϕ.A[1]) 
        + uᵣ₊₂⁻*u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻ + g*A₂) )
    # Potential energy Fᵥ
    energy += γ^4*ϕ.u⁺^2*u⁻^2*(1+ν*cos(2*(ϕ.θ⁺-ϕ.θ⁻)))
    # Andreev Bashkin terms
    energy += γ^2*(ν+1)*(u⁻*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ + g*A₂) 
        + ϕ.u⁺*uᵣ₊₂⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ + g*A₂) 
        - u⁻*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ + g*ϕ.A[1]) 
        - ϕ.u⁺*uᵣ₊₁⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ + g*ϕ.A[1]))
    # Mixed gradient terms
    energy += γ^2*(ν-1)*(ϕᵣ₊₂.u⁺*uᵣ₊₁⁻*sin(ϕᵣ₊₁.θ⁻ - ϕᵣ₊₂.θ⁺ + g*(ϕ.A[1]-A₂)) 
        - uᵣ₊₂⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕᵣ₊₂.θ⁻ + g*(ϕ.A[1] - A₂)) 
        + 2*ϕ.u⁺*u⁻*sin(ϕ.θ⁻-ϕ.θ⁺) 
        -  u⁻*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁻ + g*A₂) 
        + ϕ.u⁺*uᵣ₊₂⁻*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁺ + g*A₂) 
        + u⁻*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁻ + g*ϕ.A[1]) 
        - ϕ.u⁺*uᵣ₊₁⁻*sin(ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ + g*ϕ.A[1]))
    energy
end

# -----------------------------------------------------------------------------------------------------------
# Loops over all positions of the lattice of a state and calculates the total energy from the
# Higgs-field terms using the function fᵣ() + the energy from the gauge field.
function E(ψ::State)
    γ = ψ.γ
    g = ψ.g
    energy = 0.0
    ν = ψ.ν
    f = ψ.f
    Lx = size(ψ.lattice,2)
    Ly = size(ψ.lattice,1)
    
    # Contribution from upper right corner
    A⁰ = (Lx-1)*2π*f
    ϕ = ψ.lattice[1,Lx]    # Lattice site at upper right corner
    ϕᵣ₊₁ = ψ.lattice[1,1]    # Nearest neighbor at r+x is upper left corner
    ϕᵣ₊₂ = ψ.lattice[Ly,Lx]  # Nearest neighbor at r+y is lower right corner
    energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,γ,g,ν,f)              # Higgs terms
    energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2 # Maxwell term
    
    # Contribution from right boundary paralell to y-axis
    # except for the upper right corneϕ.
    for y=2:Ly
        ϕ = ψ.lattice[y,Lx]
        ϕᵣ₊₁ = ψ.lattice[y,1]
        ϕᵣ₊₂ = ψ.lattice[y-1,Lx]
        A⁰ = (Lx-1)*2*π*f
        energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,γ,g,ν,f)              # Higgs terms
        energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2
    end
    
    # Contribution from the bulk of lattice sites and upper boundary
    for x=1:(Lx-1)
        A⁰ = (x-1)*2*π*f        # Constant vector potential.
        # Constribution from upper boundary except upper right corner
        ϕ = ψ.lattice[1,x]
        ϕᵣ₊₁ = ψ.lattice[1,x+1]
        ϕᵣ₊₂ = ψ.lattice[Ly,x]
        energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,γ,g,ν,f)              # Higgs terms
        energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2
        
        # Contribution from the rest of the bulk.
        for y=2:Ly
            ϕ = ψ.lattice[y,x]          # Lattice site at position r
            ϕᵣ₊₁ = ψ.lattice[y,x+1]       # Nearest neighbor at r+x
            ϕᵣ₊₂ = ψ.lattice[y-1,x]       # Nearest neighbor at r+y
            
            energy += fᵣ(ϕ,ϕᵣ₊₁,ϕᵣ₊₂,A⁰,γ,g,ν,f)              # Higgs terms
            energy += (ϕ.A[1] + ϕᵣ₊₁.A[2]-ϕᵣ₊₂.A[1]-ϕ.A[2])^2
        end
    end
    
    energy
end

# -----------------------------------------------------------------------------------------------------------
# Find the energy difference between two states; one that has ϕ′ in position r with ϕᵣ... as neighbors,
# and one that has ϕ in position r. the position along the x-axis is needed for the constant Gauge field.
function ΔE(ψ::State, ϕ′::LatticeSite, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, 
        ϕᵣ₋₁::LatticeSite, ϕᵣ₋₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, ϕᵣ₋₂₊₁::LatticeSite, x::Int64)
    const g = ψ.g
    const δE::Float64 = 0.0
    
    # Calculate u⁻ for these lattice sites
    const u′⁻ = √(1-ϕ′.u⁺^2)
    const u⁻ = √(1-ϕ.u⁺^2)
    const u⁻ᵣ₊₁ = √(1-ϕᵣ₊₁.u⁺^2)
    const u⁻ᵣ₊₂ = √(1-ϕᵣ₊₂.u⁺^2)
    const u⁻ᵣ₋₁ = √(1-ϕᵣ₋₁.u⁺^2)
    const u⁻ᵣ₋₂ = √(1-ϕᵣ₋₂.u⁺^2)
    const u⁻ᵣ₋₁₊₂ = √(1-ϕᵣ₋₁₊₂.u⁺^2)
    const u⁻ᵣ₋₂₊₁ = √(1-ϕᵣ₋₂₊₁.u⁺^2)
    
    # Calculate constant link variables
    const A⁰ = 2π*ψ.f*(x-1)
    #const A⁰₊ = 2π*ψ.f*x
    const A⁰₋ = 2π*ψ.f*(x-2)
    
    # Normal kinetic terms
    δE += -2*ψ.γ^2*(ϕᵣ₊₁.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ′.θ⁺ + g*ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ + g*ϕ.A[1]))
        + ϕᵣ₋₁.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1]))
        + u⁻ᵣ₊₁*(u′⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ′.θ⁻ + g*ϕ′.A[1]) - u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ + g*ϕ.A[1]))
        + u⁻ᵣ₋₁*(u′⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1]) - u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1]))
        + ϕᵣ₊₂.u⁺*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ′.θ⁺ + g*(ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ + g*(ϕ.A[2]+A⁰)))
        + ϕᵣ₋₂.u⁺*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰)))
        + u⁻ᵣ₊₂*(u′⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ′.θ⁻ + g*(ϕ′.A[2]+A⁰)) - u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ + g*(ϕ.A[2]+A⁰)))
        + u⁻ᵣ₋₂*(u′⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰)) - u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰))))
    
    # Potential energy terms
    δE += ψ.γ^4*((ϕ′.u⁺*u′⁻)^2*(1+ψ.ν*cos(2*(ϕ′.θ⁺ - ϕ′.θ⁻))) - (ϕ.u⁺*u⁻)^2*(1+ψ.ν*cos(2*(ϕ.θ⁺ - ϕ.θ⁻))))
    
    # Andreev-Bashkin terms
    δE += ψ.γ^2*(ψ.ν+1)*(ϕᵣ₊₂.u⁺*(u′⁻*cos(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻+g*(ϕ′.A[2]+A⁰)) - u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻+g*(ϕ.A[2]+A⁰))) 
        - ϕᵣ₊₁.u⁺*(u′⁻*cos(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻+g*ϕ′.A[1]) - u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻+g*ϕ.A[1])) 
        + u⁻ᵣ₋₂*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰))) 
        - u⁻ᵣ₋₁*(ϕ′.u⁺*cos(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1]) - ϕ.u⁺*cos(ϕ.θ⁺-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1]))
        + u⁻ᵣ₊₂*(ϕ′.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺+g*(ϕ′.A[2]+A⁰)) - ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺+g*(ϕ.A[2]+A⁰))) 
        - u⁻ᵣ₊₁*(ϕ′.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺+g*ϕ′.A[1]) - ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺+g*ϕ.A[1])) 
        + ϕᵣ₋₂.u⁺*(u′⁻*cos(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰)) - u⁻*cos(ϕ.θ⁻-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰))) 
        - ϕᵣ₋₁.u⁺*(u′⁻*cos(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1]) - u⁻*cos(ϕ.θ⁻-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1])))
    
    # Mixed gradient terms
    δE += ψ.γ^2*(ψ.ν-1)*(ϕᵣ₊₂.u⁺*u⁻ᵣ₊₁*(sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ + g*(ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
            - sin(ϕᵣ₊₁.θ⁻-ϕᵣ₊₂.θ⁺ + g*(ϕ.A[1] - (ϕ.A[2]+A⁰)))) 
        -ϕᵣ₊₂.u⁺*(u′⁻*sin(ϕᵣ₊₂.θ⁺-ϕ′.θ⁻+g*(ϕ′.A[2]+A⁰)) - u⁻*sin(ϕᵣ₊₂.θ⁺-ϕ.θ⁻+g*(ϕ.A[2]+A⁰))) 
        +ϕᵣ₊₁.u⁺*(u′⁻*sin(ϕᵣ₊₁.θ⁺-ϕ′.θ⁻+g*ϕ′.A[1]) - u⁻*sin(ϕᵣ₊₁.θ⁺-ϕ.θ⁻+g*ϕ.A[1])) 
        +u⁻ᵣ₋₂₊₁*(ϕ′.u⁺*sin(ϕᵣ₋₂₊₁.θ⁻-ϕ′.θ⁺+g*(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰))) 
            - ϕ.u⁺*sin(ϕᵣ₋₂₊₁.θ⁻-ϕ.θ⁺+g*(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰)))) 
        +ϕᵣ₋₁₊₂.u⁺*(u′⁻*sin(ϕ′.θ⁻-ϕᵣ₋₁₊₂.θ⁺+g*(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋))) 
            - u⁻*sin(ϕ.θ⁻-ϕᵣ₋₁₊₂.θ⁺+g*(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋)))) 
        -u⁻ᵣ₋₂*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰)) - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₂.θ⁻+g*(ϕᵣ₋₂.A[2]+A⁰))) 
        +u⁻ᵣ₋₁*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1]) - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₁.θ⁻+g*ϕᵣ₋₁.A[1])) 
        -(u⁻ᵣ₊₂*ϕᵣ₊₁.u⁺*(sin(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻ + g*(ϕ′.A[1] - (ϕ′.A[2]+A⁰))) 
            - sin(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻ + g*(ϕ.A[1] - (ϕ.A[2]+A⁰)))) 
        -u⁻ᵣ₊₂*(ϕ′.u⁺*sin(ϕᵣ₊₂.θ⁻-ϕ′.θ⁺+g*(ϕ′.A[2]+A⁰)) - ϕ.u⁺*sin(ϕᵣ₊₂.θ⁻-ϕ.θ⁺+g*(ϕ.A[2]+A⁰))) 
        +u⁻ᵣ₊₁*(ϕ′.u⁺*sin(ϕᵣ₊₁.θ⁻-ϕ′.θ⁺+g*ϕ′.A[1]) - ϕ.u⁺*sin(ϕᵣ₊₁.θ⁻-ϕ.θ⁺+g*ϕ.A[1])) 
        +ϕᵣ₋₂₊₁.u⁺*(u′⁻*sin(ϕᵣ₋₂₊₁.θ⁺-ϕ′.θ⁻+g*(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰))) 
            - u⁻*sin(ϕᵣ₋₂₊₁.θ⁺-ϕ.θ⁻+g*(ϕᵣ₋₂.A[1] - (ϕᵣ₋₂.A[2]+A⁰)))) 
        +u⁻ᵣ₋₁₊₂*(ϕ′.u⁺*sin(ϕ′.θ⁺-ϕᵣ₋₁₊₂.θ⁻+g*(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋))) 
            - ϕ.u⁺*sin(ϕ.θ⁺-ϕᵣ₋₁₊₂.θ⁻+g*(ϕᵣ₋₁.A[1]-(ϕᵣ₋₁.A[2]+A⁰₋)))) 
        -ϕᵣ₋₂.u⁺*(u′⁻*sin(ϕ′.θ⁻-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰)) - u⁻*sin(ϕ.θ⁻-ϕᵣ₋₂.θ⁺+g*(ϕᵣ₋₂.A[2]+A⁰))) 
        +ϕᵣ₋₁.u⁺*(u′⁻*sin(ϕ′.θ⁻-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1]) - u⁻*sin(ϕ.θ⁻-ϕᵣ₋₁.θ⁺+g*ϕᵣ₋₁.A[1])))
        +2*(ϕ′.u⁺*u′⁻*sin(ϕ′.θ⁻-ϕ′.θ⁺) - ϕ.u⁺*u⁻*sin(ϕ.θ⁻-ϕ.θ⁺)))
    
    # Then calculate the Gauge field contribution
    # First contribution from current position
    δE += (ϕ′.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ′.A[2])^2 - (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])^2
    # Then from position r-x
    δE += (ϕᵣ₋₁.A[1] + ϕ′.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2 - (ϕᵣ₋₁.A[1] + ϕ.A[2] - ϕᵣ₋₁₊₂.A[1] - ϕᵣ₋₁.A[2])^2
    # Then from position r-y
    δE += (ϕᵣ₋₂.A[1] + ϕᵣ₋₂₊₁.A[2] - ϕ′.A[1] - ϕᵣ₋₂.A[2])^2 - (ϕᵣ₋₂.A[1] + ϕᵣ₋₂₊₁.A[2] - ϕ.A[1] - ϕᵣ₋₂.A[2])^2
end



####################################################################################################################
#                            Monte-Carlo functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Given a lattice site ϕ, propose a new lattice site with values in intervals around the existing ones.
function proposeLocalUpdate(ϕ::LatticeSite)
    const θmax = π/2 # How far away from the original value, the proposed new value for θ should be.
    const umax = 0.4 # How far away from the original value, the proposed amplitude u should be.
    const Amax = 3
    
    # Construct new configuration at lattice site.
    return LatticeSite([ϕ.A[1]+rand(Uniform(-Amax,Amax)), ϕ.A[2]+rand(Uniform(-Amax,Amax))],
        mod(ϕ.θ⁺ + rand(Uniform(-θmax,θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-θmax,θmax)), 2π), 
        mod(ϕ.u⁺ + rand(Uniform(-umax,umax)),1))
end

# -----------------------------------------------------------------------------------------------------------
# Performes a Metropolis Hasting update on a lattice site at position pos in state ψ given an inverse temperature
# β and where ϕᵣ... gives nearest and next nearest neighbor sites. Note that pos gives [y,x] of the position of
# the lattice site in normal array notation such that [1,1] is the upper left corner.
function metropolisHastingUpdate!(ψ::State, pos::Array{Int64,1}, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite,
        ϕᵣ₋₁::LatticeSite, ϕᵣ₋₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, ϕᵣ₋₂₊₁::LatticeSite, β::Float64)
	# Save the lattice site at the targeted position in a temporary variable ϕ and use the lattice site
	# as a basis for proposing a new lattice site ϕ′. Then find the energy difference between having
	# ϕ′ or ϕ at position pos.
    const ϕ = ψ.lattice[pos...]
    const ϕ′ = proposeLocalUpdate(ϕ)
    const δE = ΔE(ψ, ϕ′, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, pos[2])
    
    # Create random number ran ∈ (0,1].
    const ran = rand()
    if ran==0
        ran=1
    end
    
    # Update state with probability min(1, e^{-β⋅δE})
    # and return the energy of final state regardless of whether it gets updated or not.
    if log(ran) <= -β*δE
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
function mcSweep!(ψ::State, β::Float64)
   
    # Find size of the lattice L
    const L::Int64 = size(ψ.lattice,1)
    
    
    # Updating upper right corner
    const ϕᵣ₊₁ = ψ.lattice[1,1]
    const ϕᵣ₊₂ = ψ.lattice[L,L]
    const ϕᵣ₋₁ = ψ.lattice[1,L-1]
    const ϕᵣ₋₂ = ψ.lattice[2,L]
    const ϕᵣ₋₁₊₂ = ψ.lattice[L,L-1]
    const ϕᵣ₋₂₊₁ = ψ.lattice[2,1]
    metropolisHastingUpdate!(ψ, [1,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    
    # Updating boundary paralell to y-axis
    # except for the upper and lower right corner.
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,1]
        ϕᵣ₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₁ = ψ.lattice[y,L-1]
        ϕᵣ₋₂ = ψ.lattice[y+1,L]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,1]
        metropolisHastingUpdate!(ψ, [y,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    end
    
    # Updating lower right corner
    ϕᵣ₊₁ = ψ.lattice[L,1]
    ϕᵣ₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₁ = ψ.lattice[L,L-1]
    ϕᵣ₋₂ = ψ.lattice[1,L]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L-1]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,1]
    metropolisHastingUpdate!(ψ, [L,L], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    
    # Updating lower x-axis boundary except for lower left and right corner
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[L,x+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,x]
        ϕᵣ₋₁ = ψ.lattice[L,x-1]
        ϕᵣ₋₂ = ψ.lattice[1,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L-1,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[1,x+1]
        metropolisHastingUpdate!(ψ, [L,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    end
    
    # Updating lower left corner
    ϕᵣ₊₁ = ψ.lattice[L,2]
    ϕᵣ₊₂ = ψ.lattice[L-1,1]
    ϕᵣ₋₁ = ψ.lattice[L,L]
    ϕᵣ₋₂ = ψ.lattice[1,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L-1,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[1,2]
    metropolisHastingUpdate!(ψ, [L,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    
    # Updating left y-axis boundary except corners
    for y=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[y,2]
        ϕᵣ₊₂ = ψ.lattice[y-1,1]
        ϕᵣ₋₁ = ψ.lattice[y,L]
        ϕᵣ₋₂ = ψ.lattice[y+1,1]
        ϕᵣ₋₁₊₂ = ψ.lattice[y-1,L]
        ϕᵣ₋₂₊₁ = ψ.lattice[y+1,2]
        metropolisHastingUpdate!(ψ, [y,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    end
    
    # Updating upper left corner
    ϕᵣ₊₁ = ψ.lattice[1,2]
    ϕᵣ₊₂ = ψ.lattice[L,1]
    ϕᵣ₋₁ = ψ.lattice[1,L]
    ϕᵣ₋₂ = ψ.lattice[2,1]
    ϕᵣ₋₁₊₂ = ψ.lattice[L,L]
    ϕᵣ₋₂₊₁ = ψ.lattice[2,2]
    metropolisHastingUpdate!(ψ, [1,1], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    
    # Updating upper x-axis boundary
    for x=2:(L-1)
        ϕᵣ₊₁ = ψ.lattice[1,x+1]
        ϕᵣ₊₂ = ψ.lattice[L,x]
        ϕᵣ₋₁ = ψ.lattice[1,x-1]
        ϕᵣ₋₂ = ψ.lattice[2,x]
        ϕᵣ₋₁₊₂ = ψ.lattice[L,x-1]
        ϕᵣ₋₂₊₁ = ψ.lattice[2,x+1]
        metropolisHastingUpdate!(ψ, [1,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
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
        metropolisHastingUpdate!(ψ, [y,x], ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₋₁, ϕᵣ₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₂₊₁, β)
    end
end

####################################################################################################################
#                            Planar Structure function
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Returns the local vorticity by preforming a plaquette sum using the gauge-invariant difference of the θ field.
function n⁺(ψ::State, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, x::Int64)
    return (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ + ψ.g*ϕ.A[1],two_pi) + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺ + ψ.g*(ϕᵣ₊₁.A[2] + two_pi*ψ.f*x), two_pi) 
        - mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺ + ψ.g*ϕᵣ₊₂.A[1], two_pi) 
        - mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ + ψ.g*(ϕ.A[2] + two_pi*ψ.f*(x-1)), two_pi))/two_pi
end
function n⁻(ψ::State, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, x::Int64)
    return (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ + ψ.g*ϕ.A[1],two_pi) + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻ + ψ.g*(ϕᵣ₊₁.A[2] + two_pi*ψ.f*x), two_pi) 
        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻ + ψ.g*ϕᵣ₊₂.A[1], two_pi) 
        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ + ψ.g*(ϕ.A[2] + two_pi*ψ.f*(x-1)), two_pi))/two_pi
end

# -----------------------------------------------------------------------------------------------------------
function structureFunctionPluss(k::Array{Int64,1}, ψ::State)
    sum = Complex(0)
    L = size(ψ.lattice, 2)
    
    # Sum over the corners
    # Upper left corner
     r = [0, L-1] # For r we assume origo is in position [L,1] of the lattice. 
                  # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
                  # Additionally we define it such that we get the usual r = [x,y] order of dimensions.
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower left corner
     r = [0,0]
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower right corner
     r = [L-1,0]
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Upper right corner
     r = [L-1,L-1]
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        r = [i-1, L-1]
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Left border (counting in -y direction)
        r = [0, L-i]
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        
        # Lower border (counting in +x direction)
        r = [i-1, 0]
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Right border (counting in -y direction)
        r = [L-1, L-i]
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    
    return abs2(sum)
end

# -----------------------------------------------------------------------------------------------------------
# This version should do the same as above, but determines the neighbors dynamically at a slight cost
# to performance.
function structureFunctionPlussDyn(k::Array{Int64,1}, ψ::State)
    L = size(ψ.lattice, 1)
    sum = Complex(0.0)
    
    # Sum over entire lattice and determine nearest neighbors dynamically.
    for h_pos = 1:L
        for v_pos = 1:L
            
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos, mod(h_pos, L)+1]
            ϕᵣ₊₂ = ψ.lattice[mod(v_pos-2,L)+1, h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[mod(v_pos-2,L)+1, mod(h_pos, L)+1]
            sum += n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    return abs2(sum)
end



####################################################################################################################
#                            Misc. help functions
#
####################################################################################################################

# Helping function for calculating position vector in [x,y] format from a lattice position
# assuming origo is in lower left corner of the square lattice. L is the size of the lattice along one dimension.
function getVectorPosition(L::Int64, pos::Array{Int64,1})
    v_pos = pos[1]
    h_pos = pos[2]
    return [h_pos-1, L-v_pos]
end
