using Base.Test
using Distributions
using Plots
gr()

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
end
#
#   Utility functions for the types
#
import Base.copy
function copy(ϕ::LatticeSite)
    LatticeSite([ϕ.A[1],ϕ.A[2]],ϕ.θ⁺,ϕ.θ⁻,ϕ.u⁺,ϕ.u⁻)
end
function copy(c::SystConstants)
    SystConstants(c.L, c.γ, c.g⁻², c.ν, c.f, c.β)
end
function copy(ψ::State)
    Lx = size(ψ.lattice,2)
    Ly = size(ψ.lattice,1)
    lattice = [LatticeSite([ψ.lattice[y,x].A[1],ψ.lattice[y,x].A[2]],ψ.lattice[y,x].θ⁺,ψ.lattice[y,x].θ⁻,
            ψ.lattice[y,x].u⁺, ψ.lattice[y,x].u⁻) for y = 1:Ly, x=1:Lx]
    consts = copy(ψ.consts)
    State(lattice, consts)
end
function copy(sim::Controls)
    return Controls(sim.θmax, sim.umax, sim.Amax)
end

# -------------------------------------------------------------------------------------------------
# LatticeSite outer constructor
# Initializes a lattice site that has radom values
function LatticeSite()
    A_max = 3.0
    u⁺ = rand()
    LatticeSite([rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max, A_max))], two_pi*rand(), two_pi*rand(), 
        u⁺, √(1-u⁺^2))
end


####################################################################################################################
#                            Functions for ::State
#
####################################################################################################################

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
    f = 1.0/N    # Magnetic filling fraction
    β = 1/40   # Inverse temperature
    consts = SystConstants(N, γ, 1/g^2, ν, f, β)
    
    # Construct ordered state 
    if choice == 1
        
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([0,0],0,0,0,1) for y=1:N, x=1:N]
        ψ = State(lattice, consts)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
						   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), 1) for y=1:N, x=1:N]
        for y=1:N, x=1:N
            lattice[y,x].u⁻ = √(1-lattice[y,x].u⁺^2)
        end
        ψ = State(lattice, consts)
        
    # We only have choices 1 and 2 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    ψ
end
# Same as above, but now inserting the constants as SystConstants object
function State(choice::Int64, consts::SystConstants)
    N = consts.L
    N <= 1 && throw(DomainError())
    # Construct ordered state 
    if choice == 1
        
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([0,0],0,0,0,1) for y=1:N, x=1:N]
        ψ = State(lattice, consts)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
						   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), 1) for y=1:N, x=1:N]
        for y=1:N, x=1:N
            lattice[y,x].u⁻ = √(1-lattice[y,x].u⁺^2)
        end
        ψ = State(lattice, consts)
        
    # We only have choices 1 and 2 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    ψ
end

# -------------------------------------------------------------------------------------------------
# Checks if any of the lattice sites in the lattices of two states are equal and returns false in
# this case, returns true otherwise.
function latticeCompletelyDifferent(ψ₁::State, ψ₂::State)
    ψ₁.consts.L == ψ₂.consts.L || throw(DomainError())
    isCompletelyDifferent = true
    for x=1:ψ₂.consts.L, y=1:ψ₂.consts.L
        isCompletelyDifferent = isCompletelyDifferent && (ψ₁.lattice[y,x].A[1] != ψ₂.lattice[y,x].A[1]
            && ψ₁.lattice[y,x].A[2] != ψ₂.lattice[y,x].A[2] && ψ₁.lattice[y,x].θ⁺ != ψ₂.lattice[y,x].θ⁺ 
            && ψ₁.lattice[y,x].θ⁻ != ψ₂.lattice[y,x].θ⁻ && ψ₁.lattice[y,x].u⁺ != ψ₂.lattice[y,x].u⁺)
    end
    return isCompletelyDifferent
end

# -------------------------------------------------------------------------------------------------
# Checks if the state is a physically possible state. Note that we have resitricted to only having
# LxL lattices.
function checkState(ψ::State)
    N = ψ.consts.L
    @test size(ψ.lattice, 1) == size(ψ.lattice, 2) == N
    @test N>2
    for x=1:N, y=1:N
        @test ψ.lattice[y,x].θ⁺ < 2π && ψ.lattice[y,x].θ⁺ >= 0
        @test ψ.lattice[y,x].θ⁻ < 2π && ψ.lattice[y,x].θ⁻ >= 0
        @test ψ.lattice[y,x].u⁺ <= 1 && ψ.lattice[y,x].u⁺ >= 0
        @test isapprox(ψ.lattice[y,x].u⁻^2+ψ.lattice[y,x].u⁺^2, 1.0, atol=0, rtol=1e-13)
    end
    @test typeof(ψ.consts.γ) == typeof(ψ.consts.g⁻²) == typeof(ψ.consts.ν) == typeof(ψ.consts.f) == Float64
	@test isapprox(ψ.consts.L*abs(ψ.consts.f) % 1, 0.0, atol=0, rtol=1e-13)
end


####################################################################################################################
#                            Functions for ::Controls
#
####################################################################################################################

# -------------------------------------------------------------------------------------------------
function setValues!(sim_target::Controls, sim_source::Controls)
    sim_target.θmax = sim_source.θmax
    sim_target.umax = sim_source.umax
    sim_target.Amax = sim_source.Amax
    return
end


