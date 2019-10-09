# Current intention:
# Run test run using ClusterManagers to see if we can increase the CPU efficiancy

# Change list: erland <-> fram
# out_path, staging_path, src_path, split, readline comment, target temp

using Distributed
using ClusterManagers

nw = 32
time = "00:29:00"
addprocs(SlurmManager(nw), t=time)

println("After connecting workers:")
println("We have $(nprocs()) processes where workers are")
println(workers())

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

#@everywhere src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
@everywhere src_path = "/cluster/home/fredrkro/mc/Source/Grid/"
#out_path = "/home/nicolai/mc/Scripts/"
out_path = "/cluster/home/fredrkro/mc/Data/ExperimentalRuns/"
#staging_path = ""
staging_path = "/cluster/home/fredrkro/Staging/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")
include(src_path*"/test_functions.jl")

using JLD

# Enter data directory for structure function
fixRC()
# We run a simulation with parameter ν supplied by script
if length(ARGS) != 1
    println("ERROR: Need g supplied as script argument.")
    exit(-1)
end
ν = parse(Float64, ARGS[1])         # Gauge coupling
g = 0.3#[-0.5, 0.0, 0.3, 0.5, 0.7]         # Fermi surface anisotropy
κ₅ = 1.0

# TODO: We should probably benchmark vortexSnapshot and see that it is faster than mcSweeps
# Other parameters
M_est = 2^7  # Number of MC-steps to do at T_start before cooldown.
M_col = 2^7 # Number of MC-steps to use for cooling down the systems
N_steps = 2^4   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
M_th = 2^9  # Number of MC-steps to do for thermalization.
M = 2^7     # Number of measurements
M_amp = 15   # Number of these measurements that will be amplitude measurements, i.e. we need M >= M_amp
Δt = 2^6      # Number of MC-steps between each measurement
# L is assumed to be even.
L = 2*32     # System length
L₁ = L
L₂ = L
L₃ = L
split = (3,5,2)
N = L₁*L₂*L₃
# Specify extended landau gauge using n,m s.t. constant Aᵢ(Lⱼ+1) = Aᵢ(0) mod 2π
# n | L₁ and m | L₂
n = -1; m = 0
f = n/L₁ - m/L₂
T = 1.2
T_start = 2*4.0+0.1   # Start-temperature of states when thermalizing. 


# Print parameters
println("\nStarted script at: $(Dates.format(now(), "HH:MM"))")
println("Target temp: $(T)")
println("f = $(f)")
println("n = $(n)")
println("m = $(m)")
println("L = $(L)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
println("g = $(g)")
if f < 0
    println("Number of anti-vortices: $(-f*L₁*L₂)")
else
    println("Number of vortices: $(f*L₁*L₂)")
end
flush(stdout)

# Checking available workers
workers_pr_state = Int(split[1]*split[2]*split[3])
needed_workers = workers_pr_state
if nprocs()-1 < needed_workers
    println("ERROR: Only $(nprocs()-1) workers have been allocated. $(needed_workers) are needed to split
the state into $(split[1]*split[2]*split[3]) subcuboids.")
    exit()
end

# Running the state for a while
#
println("Using $(split[1]*split[2]*split[3]) workers")
syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,n,m)
cub = Cuboid(6, syst, split, 1/T_start; u⁺=1/√(2), u⁻=1/√(2), pid_start=2)

t_th = @elapsed nMCS!(cub, M*Δt)
t_MCS1 = t_th/(M*Δt)
println("Thermalization used $(round(t_th/3600; digits=1)) h, i.e. $(round(t_MCS1/(N)*1e6; digits=2)) μs pr. lattice site.")



# Checking increasing even splitting in all dir
################################################################################################
#s_list = [(1,1,1), (2,2,2), (4,2,1), (2,2,4), (2,4,4), (4,2,4), (4,4,2)]
#s_list = [(1,1,1), (4,2,4), (4,4,2)]
#for splits in s_list
#    println("Using $(splits[1]*splits[2]*splits[3]) workers")
#    syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,n,m)
#    cub = Cuboid(6, syst, splits, 1/T_start; u⁺=1/√(2), u⁻=1/√(2), pid_start=2)
#    benchmarkMcSweep!(cub; max_sec = 200)
#end



for i in workers()
    rmprocs(i)
end
