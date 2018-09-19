# Constants
const two_pi = 2π
struct SystConstants
    L::Int64      # System size in one direction
    γ::Float64    # Order parameter amplitude
    g⁻²::Float64  # 1/g² for gauge coupling g
    ν::Float64    # Anisotropy constant
    f::Float64    # Magnetic filling fraction
    β::Float64    # Simulation inverse temperature
end
mutable struct Controls
    θmax::Float64
    umax::Float64
    Amax::Float64
end

mutable struct LatticeSite
    A::Array{Float64,1}  # Fluctuating vector potential
    θ⁺::Float64 # Phase of the + component
    θ⁻::Float64 # Phase of the - component
    u⁺::Float64 # Amplitude of + component
    u⁻::Float64 # Amplitude of - component (should always be √u⁺)
end
mutable struct State
    lattice::Array{LatticeSite,2}  # Numerical lattice
    consts::SystConstants          # Collection of all constants for the state.
	nb::Array{Neighbors,2}
	nnb::Array{NextNeighbors,2}
	nnnb::Array{NNNeighbors,2}
end
struct Neighbors
	ϕᵣ₊₁::LatticeSite
	ϕᵣ₋₁::LatticeSite
	ϕᵣ₊₂::LatticeSite
	ϕᵣ₋₂::LatticeSite
end
struct NextNeighbors
	ϕᵣ₊₁₊₂::LatticeSite
	ϕᵣ₊₁₋₂::LatticeSite
	ϕᵣ₋₁₊₂::LatticeSite
	ϕᵣ₋₁₋₂::LatticeSite
end
struct NNNeighbors
	ϕᵣ₊₁₁::LatticeSite
	ϕᵣ₋₁₁::LatticeSite
	ϕᵣ₊₂₂::LatticeSite
	ϕᵣ₋₂₂::LatticeSite
end
