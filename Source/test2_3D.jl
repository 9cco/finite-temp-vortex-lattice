include("ChiralMC.jl")
using ChiralMC

#include("functions_mc.jl")
#include("functions_symmetric_energy.jl")

using Base.Test
using Distributions


L = 20
L₃ = 20
gm = 1.0
g = 0.3
ν = 0.2
κ₅ = 1.0
f = 0.03
syst = SystConstants(L, L₃, gm, g, ν, κ₅, f, 1.2)
for v=1:L, h=1:L, z=1:L₃
    ψ1 = State(1,syst)
    #Test changing one site
    ϕ = LatticeSite()

    δE = ΔE(ϕ, ψ1.lattice[v,h,z], ψ1.nb[v,h,z], ψ1.nnb[v,h,z], ψ1.nnnb[v,h,z], h, ψ1.consts)
    E0 = E(ψ1)

    set!(ψ1.lattice[v,h,z], ϕ)
    E1 = E(ψ1)
    difference = δE - (E1-E0)
    if difference > 1e-11
        println("Position: ($(v),$(h),$(z))")
        println("ΔE = $(δE), E1-E0 = $(E1-E0), E1 = $(E1), E0 = $(E0)")
        
        println("difference: $(δE - (E1-E0))")
    end
#    println(@test isapprox(δE, E1-E0, atol = 0, rtol = 1e-13))
end


