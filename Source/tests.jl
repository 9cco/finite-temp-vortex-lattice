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
println(@test size(ψ.lattice,1) == size(ψ.lattice,2) && size(ψ.lattice,2) == N)
println("Checking that ψ is a valid State")
function checkState(ψ::State)
    for x=1:N, y=1:N
        @test ψ.lattice[y,x].θ⁺ < 2π && ψ.lattice[y,x].θ⁺ >= 0
        @test ψ.lattice[y,x].θ⁻ < 2π && ψ.lattice[y,x].θ⁻ >= 0
        @test ψ.lattice[y,x].u⁺ <= 1 && ψ.lattice[y,x].u⁺ >= 0
    end
    @test typeof(ψ.γ) == typeof(ψ.g) == typeof(ψ.ν) == typeof(ψ.f) == Float64
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
