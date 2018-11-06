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

#end

#for v_pos=1:L, h_pos=1:L, z_pos=1:L₃
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
    nnnbl


#########################################################################################
#       Testing energy functions for 3D case
#
#########################################################################################

println("\nTesting energy functions\n")
#First test the kinetic and potential energies for a state with no gauge field
# fluctuations and phases set to zero

L = 100
L₃ = 60
f = 1/L*(rand()+1)
gm = 1.0
g = 0.3
ν = rand()
κ₅ = 1.0

syst = SystConstants(L, L₃, gm, g, ν, κ₅, f, 0.0)
ψ = State(1, syst)




