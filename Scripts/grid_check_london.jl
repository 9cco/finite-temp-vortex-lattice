# Current intention:
# Run grid-parallelized code on a one-component london system for a series of different temperatures to see that the dual
# stiffness is as expected and has a kink at T ≈ 1.9

using Distributed
@everywhere using Distributions
@everywhere using Test
using BenchmarkTools
using Dates

@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()

src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere src_path = "../Source/Grid/"
@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")

using Plots
pyplot()
using JLD

# Enter data directory for structure function
fixRC()
# We run a simulation with the parameters
g = √(1/10)    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
M_col = 2^9#2^17 # Number of MC-steps to use for cooling down the systems
M_th = 0#2^13  # Number of MC-steps to do for thermalization.
M = 2^2#2^13    # Number of measurements
Δt = 19    # Number of MC-steps between each measurement
# L is assumed to be even.
L = 24     # System length
L₁ = L
L₂ = L
L₃ = L
N = L₁*L₂*L₃
# Make geometric progression of temperatures between T₁ and Tₘ
T₁ = 0.7
Tₘ = 4.0
N_T = 3*2^0
N_steps = 2^3   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
T_start = 2*4.0+0.1   # Start temperature of states when thermalizing. Must be higher than maximum(temps)
R = (Tₘ/T₁)^(1/(N_T-1))
temps = 2*[0.35, 0.4, 1.66]#[T₁*R^(k-1) for k = 1:N_T]#
println("Temps.: $(temps)")
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 0.0/L
println("f = $(f)")
println("L = $(L)")
println("g = $(g)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
mkcd("one_comp_london_L=$(L)_T=$(round(temps[1]; digits=2))-$(round(temps[end]; digits=2))_g=$(round(g; digits=3))_fL=$(round(f*L; digits=1))_grid_test")

# Make ab inito un-correlated phases state
init_syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,f)
cubs = [Cuboid(1, init_syst, (1,1,3), 1/T_start; u⁺=0.0) for T in temps]


# Time-estimation
################################################################################################

M_est = 2^5
# Do M_est MC-sweeps on first
t_meas = @elapsed for i = 1:M_est
    for j = 1:Δt
        for cub in cubs
            mcSweepEnUp!(cub)
        end
    end

    all_res = [dualStiffnesses(cub) for cub in cubs]
    En_res = [energy(cub) for cub in cubs]
end
t_meas = t_meas/M_est
t_MCS = t_meas/(Δt*N_T)
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS")
println("and $(round(t_MCS/N*1e6; digits=1)) μs on average pr. lattice site.")

# Estimating ETC
println("Measurements and thermalization will do $(Δt*M + M_th + M_col) MCS for each of the $(N_T)
temperatures, which has an ETC of $(round((Δt*M + M_th + M_col)*N_T*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")
flush(stdout)

user_input = readline(stdin)
if user_input == "n"
    exit()
end


# Cooldown
################################################################################################
println("\nStarted cooldown at: $(Dates.format(now(), "HH:MM"))")

# Do M_th MCS in N_steps steps
M_pr_step = ceil(Int64, M_col/N_steps)
M_col = M_pr_step*N_steps
# Each step is associated with a different set of temperatures given by each row in the matrix
temp_mt = genGeometricTemperatureSteps(T_start, temps, N_steps)
E_matrix = Array{Float64}(undef, M_col+1, N_T);
for (i, cub) = enumerate(cubs); E_matrix[1,i] = energy(cub); end
E_therm = Array{Float64, 2}(undef, M_th, N_T);

println("Cooling down from T=$(T_start) using $(N_steps) steps with $(M_pr_step) MCS pr. step.")
flush(stdout)

t_col = @elapsed for step = 1:N_steps
    # First set the temperatures associated with this step to each system.
    β_step = [1/T for T in temp_mt[step, :]]
    for (i, cub) = enumerate(cubs); setTemp!(cub, temp_mt[step, i]); end

    # Do the MCS associated with this step
    for i = 1:M_pr_step
        for (j, cub) = enumerate(cubs)
            δE, _ = mcSweepEnUp!(cub)
            E_matrix[(step-1)*M_pr_step+i+1, j] = E_matrix[(step-1)*M_pr_step+i, j] + δE
        end
    end
    # In the end of a step we correct energy rounding-error
    for (j, cub) = enumerate(cubs); E_matrix[step*M_pr_step+1, j] = energy(cub); end
end

E_matrix = E_matrix./N;

int = 1:M_col+1#floor(Int64, M_th/2)
therm_plt = plot(collect(int).+M_est, [E_matrix[int, i] for i = 1:N_T]; 
                 label=reshape(["T = $(round(T; digits=2))" for T in temps], (1, N_T)),
                 xaxis="MCS", yaxis="Energy pr. site", title="Cooldown energy from T=$(T_start)")
savefig(therm_plt, "cooldown energies.pdf")



# Adjusting simulation constants
################################################################################################

println("Adjusting simulation constants.")
@time for cub in cubs; tune!(cub); end
println("Controls after adjustment:")
printControls(cubs)


