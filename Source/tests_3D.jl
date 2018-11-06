include("ChiralMC.jl")
using ChiralMC

include("functions_mc.jl")
include("functions_observables.jl")
include("functions_symmetric_energy.jl")

using Base.Test
using Distributions

#########################################################################################
#       Testing neighbor functions
#
#########################################################################################

println("Testing neighbor function")

L = 5
L₃ = 4
f = 1/L*(rand()+1)
gm = 1.0
g = 0.3
ν = rand()
κ₅ = 1.0

syst = SystConstants(L, L₃, gm, g, ν, κ₅, f, 0.0)
ψ = State(1, syst)

println("Testing that all neighbors are asigned")
for h_pos=1:L, v_pos = 1:L, z_pos = 1:L₃
    @test isassigned(ψ.nb, v_pos, h_pos, z_pos)
    @test isassigned(ψ.nnb, v_pos, h_pos, z_pos)
    @test isassigned(ψ.nnnb, v_pos, h_pos, z_pos)
end

println("\nTesting lattice neighbors")
#Use phases and + comp as dummy indices for the positions
lattice = [LatticeSite([0,0],v_pos,h_pos,z_pos,1) for v_pos=1:L, h_pos=1:L, z_pos=1:L₃]
nbl = latticeNeighbors(lattice,L,L₃)
nnbl = latticeNextNeighbors(lattice,L,L₃)
nnnbl = latticeNNNeighbors(lattice,L,L₃)


#TODO: Sjekke dette bedre i morgen (sjekket og funker)
for v_pos=1:L, h_pos=1:L, z_pos=1:L₃
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₁.u⁺ == Float64(z_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₁.u⁺ == Float64(z_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₂.u⁺ == Float64(z_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₂.u⁺ == Float64(z_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₃.u⁺ == Float64(mod(z_pos-2,L₃)+1)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₃.u⁺ == Float64(mod(z_pos,L₃)+1)

    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₃.θ⁺ == Float64(v_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₃.θ⁺ == Float64(v_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₃.θ⁻ == Float64(h_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₃.θ⁻ == Float64(h_pos)
end

println("Testing latticeNeighbors")
for v_pos=1:L, h_pos=1:L, z_pos = 1:L₃
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₁.θ⁺ == Float64(v_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₁.θ⁻ == Float64(mod(h_pos,L)+1)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₁.θ⁺ == Float64(v_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₁.θ⁻ == Float64(mod(h_pos-2,L)+1)
    
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₂.θ⁺ == Float64(mod(v_pos-2,L)+1)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₊₂.θ⁻ == Float64(h_pos)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₂.θ⁺ == Float64(mod(v_pos,L)+1)
    @test nbl[v_pos, h_pos, z_pos].ϕᵣ₋₂.θ⁻ == Float64(h_pos)
end #Fra gammel kode, disse funker

#Random stikkprøver
@test nbl[L,L,L₃].ϕᵣ₋₃.u⁺ == 1
@test nbl[L,L,L₃].ϕᵣ₋₂.θ⁺ == 1
@test nbl[L,L,L₃].ϕᵣ₊₂.θ⁺ == L-1
@test nbl[L,L,L₃].ϕᵣ₊₃.u⁺ == L₃-1
@test nbl[L,L,L₃].ϕᵣ₋₁.θ⁻ == L-1
@test nbl[L,L,L₃].ϕᵣ₊₁.θ⁻ == 1

#TODO: Test for NNeighbors i z-retning, må implenteres før de kan testes
#println("\nTesting LaticeNextNeighbors")
#for v_pos=1:L, h_pos=1:L, z_pos=1:L₃
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₊₃.u⁺ == Float64(mod(z_pos-2,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₋₃.u⁺ == Float64(mod(z_pos,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₊₃.u⁺ == Float64(mod(z_pos-2,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₋₃.u⁺ == Float64(mod(z_pos,L₃)+1)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₊₃.u⁺ == Float64(mod(z_pos-2,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₋₃.u⁺ == Float64(mod(z_pos,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₊₃.u⁺ == Float64(mod(z_pos-2,L₃)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₋₃.u⁺ == Float64(mod(z_pos,L₃)+1)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₊₂.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₋₂.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₊₂.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₋₂.u⁺ == Float64(z_pos)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₊₁.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₋₁.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₊₁.u⁺ == Float64(z_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₋₁.u⁺ == Float64(z_pos)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₊₃.θ⁺ == Float64(v_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₋₃.θ⁺ == Float64(v_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₊₃.θ⁺ == Float64(mod(v_pos-2,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₋₃.θ⁺ == Float64(mod(v_pos-2,L)+1)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₊₃.θ⁺ == Float64(v_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₋₃.θ⁺ == Float64(v_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₊₃.θ⁺ == Float64(mod(v_pos,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₋₃.θ⁺ == Float64(mod(v_pos,L)+1)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₊₃.θ⁻ == Float64(mod(h_pos,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₋₃.θ⁻ == Float64(mod(h_pos,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₊₃.θ⁻ == Float64(h_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₋₃.θ⁻ == Float64(h_pos)

#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₊₃.θ⁻ == Float64(mod(h_pos-2,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₋₃.θ⁻ == Float64(mod(h_pos-2,L)+1)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₊₃.θ⁻ == Float64(h_pos)
#    @test nnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₋₃.θ⁻ == Float64(h_pos)
#end


println("\nTesting LatticeNNNeighbors")
for v_pos=1:L, h_pos=1:L, z_pos=1:L₃
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₊₃₃.u⁺ == Float64(mod(z_pos-3,L₃)+1)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₋₃₃.u⁺ == Float64(mod(z_pos+1,L₃)+1)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₊₃₃.θ⁺ == Float64(v_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₋₃₃.θ⁺ == Float64(v_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₊₃₃.θ⁻ == Float64(h_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₋₃₃.θ⁻ == Float64(h_pos)

    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₊₁₁.u⁺ == Float64(z_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₋₁₁.u⁺ == Float64(z_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₊₂₂.u⁺ == Float64(z_pos)
    nnnbl[v_pos, h_pos, z_pos].ϕᵣ₋₂₂.u⁺ == Float64(z_pos)
end

#########################################################################################
#       Testing energy functions for 3D case
#
#########################################################################################

println("\nTesting energy functions\n")
#First test the kinetic and potential energies for a state with no gauge field
# (no fluctuating field and no external)
# fluctuations and phases set to zero

L = 5
L₃ = 4
f = 0.0
gm = 0.3
g = 0.3
ν = 0.5
κ₅ = 1.0

syst = SystConstants(L, L₃, gm, g, ν, κ₅, f, 0.0)
ψ = State(1, syst)

energy_theoretical = L^2*L₃*(-gm^2 + gm^4/2) 
println(@test isapprox(energy_theoretical, E(ψ), atol = 0, rtol = 1e-13))

#Set the u⁺ amplitudes to one as well
for v_pos=1:L, h_pos=1:L, z_pos = 1:L₃
    ψ.lattice[v_pos,h_pos,z_pos].u⁺ = 1.0
end

energy_theoretical = L*L*L₃*gm^2*(gm^2*(3+ν)-2)
println(@test isapprox(energy_theoretical, E(ψ), atol = 0, rtol = 1e-13))

#Non-zero angles 
θ⁺ = π/2
θ⁻ = π/7
for v_pos=1:L, h_pos=1:L, z_pos=1:L₃
    ψ.lattice[v_pos,h_pos,z_pos].θ⁺ = θ⁺
    ψ.lattice[v_pos,h_pos,z_pos].θ⁻ = θ⁻
end

energy_theoretical = L*L*L₃*gm^2*(gm^2*(3+ν*cos(2*(θ⁺-θ⁻)))-2)
println(@test isapprox(energy_theoretical, E(ψ), atol = 0, rtol = 1e-13))

#Include a constant gauge field (uniform background), only u⁻ non-zero otherwise.
f = 1/L*(rand()+1)

syst = Systconstants(L,L₃,gm,g,ν,κ₅,f,0.0)
ψ = State(1,syst)



