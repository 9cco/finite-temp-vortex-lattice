include("types.jl")
include("functions.jl")

using Base.Test

# Testing State constructor.
const N = 20
const M = 100
const β = 1/300
ψ = State(N, 2)

println("Testing State constructor\n----------------------------------------------------------------")
println("Testing if lattice dimensions are correct")
println(@test size(ψ.lattice,1) == size(ψ.lattice,2) && size(ψ.lattice,2) == N && N == ψ.consts.L)
println("Checking that ψ is a valid State")
function checkState(ψ::State)
    for x=1:N, y=1:N
        @test ψ.lattice[y,x].θ⁺ < 2π && ψ.lattice[y,x].θ⁺ >= 0
        @test ψ.lattice[y,x].θ⁻ < 2π && ψ.lattice[y,x].θ⁻ >= 0
        @test ψ.lattice[y,x].u⁺ <= 1 && ψ.lattice[y,x].u⁺ >= 0
        @test isapprox(ψ.lattice[y,x].u⁻^2+ψ.lattice[y,x].u⁺^2, 1.0, atol=0, rtol=1e-13)
    end
    @test typeof(ψ.consts.γ) == typeof(ψ.consts.g) == typeof(ψ.consts.ν) == typeof(ψ.consts.f) == Float64
end
println(checkState(ψ))
# Testing fᵣ


# Testing mcSweep!
println("\nTesting mcSweep!\n----------------------------------------------------------------")
ψ_old = copy(ψ)
println("Checking if copied state is correct")
println(checkState(ψ_old))

# Do mcSweep! M times
for i = 1:M
    mcSweep!(ψ, β)
end
println("Checking if mcSweeped state is a state")
println(checkState(ψ))
println("Checking if mcSweeped state has lower energy")
println(@test E(ψ) < E(ψ_old))


########################################################################################
#                     Testing Energy functions
#
#######################################################################################
println("\nTesting Energy function E(ψ)\n----------------------------------------------------------------")
# Test whether the normal kinetic energy terms get correctly calculated when all fluctuating gauge fields
# vanish as well as all phases. In this case the expression for the free energy simplifies

L = 100
f = 1/L*(rand()+1)
γ = 1.0
ν = rand()

# First we create a state with all phases and fluctuating gauge fields 0.
ψ = State(L, 1, γ, 1.0, ν, f)

println("Checking that kinetic energy is calculated correctly")
# Then we test if the kinetic energy is correctly calculated as
# Fₖ = -2γ^2(L^2 - L*sin(πfL)cos(πf(L-1))/sin(πf))
@test isapprox(E(ψ), -2*γ^2*L*(L+sin(π*f*L)*cos(pi*f*(L-1))/sin(π*f)), atol=0, rtol=1e-13)


# Test whether kinetic energy terms, potential terms and Andreev Bashkin terms are calculated correctly when
# fluctuating gauge fields vanish as well as phases, while u⁺=u⁻=1/√2
# In this case, these energies should evaluate to
#
# F_K = -γ²[L² + L⋅χ]
# F_V = γ⁴L²(1+ν)/2⁴
# F_AB = γ²(ν+1)/2⋅(L⋅χ - L²)
#
# where χ = sum_{x=0}^{L-1} cos(2πfx)

L = 100
f = 1/L*(rand()+1)
γ = 1.0
ν = rand()

# Create a state where u⁺=u⁻=1/√2 while Aᵣᵤ = δᵤ₂⋅2πfrₓ and θ⁺=θ⁻=0
ψ = State(L, 1, γ, 1.0, ν, f)
for x = 1:L
    for y = 1:L
        ψ.lattice[y,x].u⁺ = 1/sqrt(2)
    end
end

# Calculate the theoretical energy
χ = sum([cos(two_pi*f*x) for x=0:(L-1)])
teoretisk = γ^2*L*(χ*(ν-1)-L*(ν+3) + γ^2*L*(1+ν)/4)

println("Checking that potential terms and Andreev Bashkin terms are calculated correctly")
@test isapprox(E(ψ),teoretisk,atol=0,rtol=1e-13)


# Test whether the Maxwell terms
# are calculated correctly by creating a system where u⁺=0, u⁻=1, θ⁺=θ⁻=0 and the 
# fluctuating gauge fields vanish except for a plaquette at [2,2].
# In this case, the energies that contribute are the maxwell terms and the kinetic energy which becomes
#
# F_K = -2γ²{ cos(Aᵣ,₁) + cos(Aᵣ₊₁,₂ + 2πf⋅2) + cos(Aᵣ₊₂,₁) + cos(Aᵣ,₂ + 2πf⋅1) + L²-2 + L⋅χ - cos(2πf⋅2) - cos(2πf⋅1) }
# F_A = -1/g⋅[ (Aᵣ,₁ + Aᵣ₊₁,₂ - Aᵣ₊₂,₁ - Aᵣ,₂)² + Aᵣ,₁² + Aᵣ₊₁,₂² + Aᵣ₊₂,₁² + Aᵣ,₂² ]
#
# where χ = sum_{x=0}^{L-1} cos(2πfx)

L = 100
f = 1/L*(rand()+1)
γ = 1.0
ν = rand()
g = 1.0

# Create a state where u⁺=u⁻=1/√2 while Aᵣᵤ = δᵤ₂⋅2πfrₓ and θ⁺=θ⁻=0
A_max = 3
ψ = State(L, 1, γ, g, ν, f)
# Creating random numbers at a plaquette
ψ.lattice[2,2].A[1] = rand(Uniform(-A_max,A_max))
ψ.lattice[2,3].A[2] = rand(Uniform(-A_max,A_max))
ψ.lattice[1,2].A[1] = rand(Uniform(-A_max,A_max))
ψ.lattice[2,2].A[2] = rand(Uniform(-A_max,A_max))

# Calculate the theoretical energy
χ = sum([cos(two_pi*f*x) for x=0:(L-1)])
teoretisk = -2*γ^2*( cos(ψ.lattice[2,2].A[1]) + cos(ψ.lattice[2,3].A[2] + 2π*f*2) + cos(ψ.lattice[1,2].A[1]) 
    + cos(ψ.lattice[2,2].A[2] + 2π*f) + L^2-2 + L*χ - cos(2π*f*2) - cos(2π*f) )
teoretisk += -( (ψ.lattice[2,2].A[1] + ψ.lattice[2,3].A[2] - ψ.lattice[1,2].A[1] - ψ.lattice[2,2].A[2])^2 
    + ψ.lattice[2,2].A[1]^2 + ψ.lattice[2,3].A[2]^2 + ψ.lattice[1,2].A[1]^2 + ψ.lattice[2,2].A[2]^2 )/g

println("Checking that fluctuating gauge field is handeled correctly for test case")
println(@test isapprox(E(ψ),teoretisk,atol=0,rtol=1e-13))

E_old = E(ψ)
ψ.g = g+10
println("Checking that increasing coupling constant g increases the Maxwell term energy")
println(@test E(ψ) > E_old)



println("\nTesting Energy difference function ΔE\n----------------------------------------------------------------")
# Test case. We make a random 3x3 lattice, and put two different lattice sites in the middle
# i.e. the [2,2] position. Then we use the different functions to calculate the energy
# difference associated with this change.

ψ₂ = State(3,2)
site = LatticeSite() # Get random lattice site
ψ₁ = copy(ψ₂)
ψ₁.lattice[2,2] = copy(site)

dE = ΔE(ψ₂,site,ψ₂.lattice[2,2],ψ₂.lattice[2,3],ψ₂.lattice[1,2],ψ₂.lattice[2,1],ψ₂.lattice[3,2],ψ₂.lattice[1,1]
    ,ψ₂.lattice[3,3], 2)
println("Checking that ΔE and E get same result")
@test isapprox(E(ψ₁)-E(ψ₂), dE; atol=0, rtol=1e-13)

########################################################################################
#                     Testing local vorticity
#
#######################################################################################
println("\nTesting n⁺\n----------------------------------------------------------------")
# Creating state where there is no contribution from Gauge fields
ψ = State(2, 2, 4.0, 0.0, 1.0, 0.0)
ϕ = ψ.lattice[2,1]
ϕᵣ₊₁ = ψ.lattice[2,2]
ϕᵣ₊₂ = ψ.lattice[1,1]
ϕᵣ₊₁₊₂ = ψ.lattice[1,2]
println("Checking n⁺ for difference phase configurations.")
# Setting phases manually
ϕ.θ⁺ = π/4
ϕᵣ₊₁.θ⁺ = 3π/4
ϕᵣ₊₂.θ⁺ = 3π/4
ϕᵣ₊₁₊₂.θ⁺ = 5π/4
println(@test n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 0) == 0)
ϕ.θ⁺ = π/4
ϕᵣ₊₁.θ⁺ = 3π/4
ϕᵣ₊₂.θ⁺ = 7π/4
ϕᵣ₊₁₊₂.θ⁺ = 5π/4
println(@test n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 0) == -1)
ϕ.θ⁺ = π/4
ϕᵣ₊₁.θ⁺ = -π/4
ϕᵣ₊₂.θ⁺ = -5π/4
ϕᵣ₊₁₊₂.θ⁺ = -3π/4
println(@test n⁺(ψ, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 0) == 1)


########################################################################################
#                     Testing structure function sum
#
#######################################################################################
println("\nTesting structureFunctionPluss(k, ψ)\n----------------------------------------------------------------")
# First the zero test. If we set everything to 0 do we get 0 in the end?
ψ = State(30, 1)
L = size(ψ.lattice,1)
ψ.f = 0
k = [rand(1:L)-1, rand(1:L)-1]
println(@test structureFunctionPluss(k, ψ) == 0)

# Even if we create a local vortex at [L-1,2], the sum should be zero still if k=0
ψ.lattice[L-1, 3].θ⁺ = π/2
ψ.lattice[L-2, 3].θ⁺ = π
ψ.lattice[L-2, 2].θ⁺ = 3π/2
println(@test n⁺(ψ,ψ.lattice[L-1,2],ψ.lattice[L-1,3],ψ.lattice[L-2,2],ψ.lattice[L-2,3], 2) == -1)
println(@test structureFunctionPluss([0,0], ψ) == 0)

# For a general k, we will get the expression
# Note that the order of the terms in the sum is important for getting the same floating point value.
res = 0.0
res -= exp(im*(k⋅getVectorPosition(L,[L,3])))
res += exp(im*(k⋅getVectorPosition(L,[L-2, 1])))
res += exp(im*(k⋅getVectorPosition(L,[L-2, 2])))
res -= exp(im*(k⋅getVectorPosition(L,[L-1,2])))
res = abs2(res)
println(@test structureFunctionPluss(k, ψ) == res)

# Finally we want to turn on the static part of the gauge field and make sure that this normalizes correctly
# when k = 0, given that all phases are again zero.
ψ = State(30,1)
ψ.f = 1/(2*L)
ψ.g = 1
println(@test isapprox(structureFunctionPluss([0,0], ψ),(ψ.f*L^2)^2,atol=0,rtol=1e-13))
