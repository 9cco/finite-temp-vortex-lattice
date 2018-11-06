# Constants
const two_pi = 2π
struct SystConstants
    L::Int64      # System size in one direction
    L₃::Int64     # System size in z-direction
    γ::Float64    # Order parameter amplitude
    g⁻²::Float64  # 1/g² for gauge coupling g
    ν::Float64    # Anisotropy constant
    κ₅::Float64   # z-layer coupling.
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
struct NearestNeighbors
	ϕᵣ₊₁::LatticeSite
	ϕᵣ₋₁::LatticeSite
	ϕᵣ₊₂::LatticeSite
	ϕᵣ₋₂::LatticeSite
    ϕᵣ₊₃::LatticeSite
    ϕᵣ₋₃::LatticeSite
end
struct NextNeighbors
	ϕᵣ₊₁₊₂::LatticeSite
	ϕᵣ₊₁₋₂::LatticeSite
	ϕᵣ₋₁₊₂::LatticeSite
	ϕᵣ₋₁₋₂::LatticeSite
    ϕᵣ₊₁₋₃::LatticeSite
    ϕᵣ₋₁₊₃::LatticeSite
    ϕᵣ₊₂₋₃::LatticeSite
    ϕᵣ₋₂₊₃::LatticeSite
end
struct NNNeighbors
	ϕᵣ₊₁₁::LatticeSite
	ϕᵣ₋₁₁::LatticeSite
	ϕᵣ₊₂₂::LatticeSite
	ϕᵣ₋₂₂::LatticeSite
    ϕᵣ₊₃₃::LatticeSite
    ϕᵣ₋₃₃::LatticeSite
end
struct Neighbors
    nb::NearestNeighbors
    nnb::NextNeighbors
    nnnb::NNNeighbors
end
mutable struct State
    lattice::Array{LatticeSite,3}  # Numerical lattice
    consts::SystConstants          # Collection of all constants for the state.
	nb::Array{NearestNeighbors,3}
	nnb::Array{NextNeighbors,3}
	nnnb::Array{NNNeighbors,3}
end
