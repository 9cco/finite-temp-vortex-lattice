# Current intention:
# Load the code and run some benchmarks to see if we can use more cores etc.
# We assume that the script is called with -p 31 s.t. nprocs()==32. We estimate
# the time it takes to run mcSweep!(cub) pr lattice site for an increasing number 
# of splittings of the lattice.

using Distributed
@everywhere using Distributions
using Test
using BenchmarkTools
using Dates

@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()

@everywhere src_path = "/cluster/projects/nn2819k/finite-temp-vortex-lattice/Source/Grid"
#@everywhere src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

#@everywhere include(src_path*"/observables.jl")
include(src_path*"/test_functions.jl")

fixRC()
# We run a simulation with the parameters
g = 0.3# √(1/10)    # Gauge coupling
ν = 0.3    # Anisotropy

# TODO: We should probably benchmark vortexSnapshot and see that it is faster than mcSweeps
# Other parameters
# L is assumed to be even.
T = 2*4.0+0.1 
κ₅ = 1.0

println("g = $(g)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
flush(stdout)



# Checking increasing splitting in z-direction
################################################################################################


# Checking increasing splitting in x-direction
################################################################################################

# Given a (2,2,2) splitting, check different lattice sizes
################################################################################################

# Checking increasing even splitting in all dir
################################################################################################
s_list = [(1,1,1), (1,1,7), (1,1,8), (2,2,2), (1,2,4), (2,1,4), (4,2,1)]
L₁ = L₂ = L₃ = 32
f = 1.0/L₁
for splits in s_list
    syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,f)
    cub = Cuboid(1, syst, splits, 1/T; u⁺=0.0)
    benchmarkMcSweep!(cub; max_sec = 200)
end

# Given equal shell-to-lattice-site ratio by changing L and splitting, checking for increasing splitting
################################################################################################

# Time-estimation
################################################################################################

