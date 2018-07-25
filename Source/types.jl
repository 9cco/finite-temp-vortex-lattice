type LatticeSite
    A::Array{Float64,1}  # Fluctuating vector potential
    θ⁺::Float64 # Phase of the + component
    θ⁻::Float64 # Phase of the - component
    u⁺::Float64 # Amplitude of + component
end
type State
    lattice::Array{LatticeSite,2}  # Numerical lattice
    γ::Float64    # Order parameter amplitude
    g::Float64    # Gauge coupling
    ν::Float64    # Anisotropy constant
    f::Float64    # Magnetic filling fraction
end
#
#   Utility functions for the types
#
import Base.copy
function copy(ϕ::LatticeSite)
    LatticeSite([ϕ.A[1],ϕ.A[2]],ϕ.θ⁺,ϕ.θ⁻,ϕ.u⁺)
end
function copy(ψ::State)
    Lx = size(ψ.lattice,2)
    Ly = size(ψ.lattice,1)
    lattice = [LatticeSite([ψ.lattice[y,x].A[1],ψ.lattice[y,x].A[2]],ψ.lattice[y,x].θ⁺,ψ.lattice[y,x].θ⁻,
            ψ.lattice[y,x].u⁺) for y = 1:Ly, x=1:Lx]
    State(lattice, ψ.γ, ψ.g, ψ.ν, ψ.f)
end

# -------------------------------------------------------------------------------------------------
# Outer constructor
# Initializes a state ψ that either 1: has zero as value for the the fluctuating gauge potential link variables,
# the phase and the u⁺ component (which means the u⁻=1) at each lattice site, or 2: has random values for
# these variables. The lattice in state ψ will consist of NxN lattice sites.
function State(N::Int64, choice::Int64)
    N <= 1 && throw(DomainError())
    
    # Constants
    γ = 1.0    # Order parameter amplitude
    g = 1.0    # Gauge coupling
    ν = 0.0    # Anisotropy constant
    f = 0.0    # Magnetic filling fraction
    
    # Construct ordered state 
    if choice == 1
        
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([0,0],0,0,0) for y=1:N, x=1:N]
        ψ = State(lattice, γ, g, ν, f)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
							   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand()) for y=1:N, x=1:N]
        ψ = State(lattice, γ, g, ν, f)
        
    # We only have choices 1 and 2 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    ψ
end

