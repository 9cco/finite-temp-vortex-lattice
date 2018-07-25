include("types.jl")
include("functions.jl")

using Base.Test

# Testing State constructor.
const N = 20
const M = 100
const β = 1/300
const ψ = State(N, 2)

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
