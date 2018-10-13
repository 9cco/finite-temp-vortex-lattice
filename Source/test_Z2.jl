#Quick test of function to measure Z2 order parameter

println("Loading MC-code")

@everywhere include("ChiralMC.jl")
@everywhere using ChiralMC
@everywhere include("functions_parallel.jl")
@everywhere include("functions_observables.jl")

#Simulation parameters
L=20
β=10^4
gm=1.0
g=0.3
ν=0.2
f=0
c = SystConstants(L, gm, 1/g^2, ν, f, β)
sim = Controls(π/3, 0.4, 3.0)

n_workers = 1

ψ_ref = State(2,L)
ψ_w = [State(1,L) for i=1:n_workers]



MCS, E_ref, E_w, ψ_ref, ψ_w, sim_ref, sim_w = parallelThermalisation!(ψ_ref, ψ_w, c, sim, 1000, 1.5, 0.1)



Z2_measure = Array{Complex}(L,L)
Z₂MeasureLocal!(ψ_w[1],Z2_measure)
#k = [0.0,0.0]
#Z₂Fourrier=Z₂StructureFunction(k, ψ_w[1], Z2_measure, 10, 10, L)

#println(Z2_measure)
#println(abs(mean(Z2_measure)))
#println(Z₂Fourrier)



using PyPlot
N=101
K = [-π + 2*π*i/(N-1) for i = 0:N-1]
Z₂FourrierTrans = Array{Complex}(N,N)
for i=1:N
    for j=1:N
        Z₂FourrierTrans[i,j] = Z₂StructureFunction([K[i],K[j]], ψ_w[1], Z2_measure,10,10,L)
    end
end

Z₂FtransAbs = abs.(Z₂FourrierTrans)
p1=surf(K,K,Z₂FtransAbs)
show()


#Energy = 0.0
#Z2_measure = Complex(0)
#MaxMCS = 100000
#Measure_period = 10
#for i=1:MaxMCS
#     Energy += mcSweepEn!(ψ_ref, sim_ref)
#    if i % Measure_period == 0
#        Z2_measure += Z2(ψ_ref)
#    end
#end
#Z2_measure = Z2_measure*Measure_period/MaxMCS
#println("Z2 measured: $(Z2_measure)")
#println("Absolute Z2 measured: $(abs(Z2_measure))")
#println("Energy change during MCS measure steps: $(Energy)")
    
