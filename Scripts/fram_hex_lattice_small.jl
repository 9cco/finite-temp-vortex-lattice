# Current intention:
# Running simulation on Fram. Try a small version of large system to check that everything is a-ok.

using Distributed
@everywhere using Distributions
using Test
#using BenchmarkTools
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
out_path = "/cluster/home/fredrkro/mc/Data/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")

#using Plots
#pyplot()
using JLD

# Enter data directory for structure function
fixRC()
# We run a simulation with the parameters
g = 0.1# √(1/10)    # Gauge coupling
ν = 0.3    # Anisotropy

# TODO: We should probably benchmark vortexSnapshot and see that it is faster than mcSweeps
# Other parameters
M_est = 2^5 # Number of MC-steps to do at T_start before cooldown.
M_col = 2^13 # Number of MC-steps to use for cooling down the systems
M_th = 2^12#2^13  # Number of MC-steps to do for thermalization.
M = 2^12#2^13    # Number of measurements
Δt = 2^2#14    # Number of MC-steps between each measurement
# L is assumed to be even.
L = 32     # System length
L₁ = L
L₂ = L
L₃ = L
split = (2,2,2)
N = L₁*L₂*L₃
# Make geometric progression of temperatures between T₁ and Tₘ
T₁ = 1.7
Tₘ = 2.0
N_T = 1#3*2^0
N_steps = 2^6   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
T_start = 2*4.0+0.1   # Start temperature of states when thermalizing. Must be higher than maximum(temps)
R = (Tₘ/T₁)^(1/(N_T-1))
temps = [0.8]#2*[0.35, 0.4, 1.66]#[T₁*R^(k-1) for k = 1:N_T]#
println("Temps.: $(temps)")
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 3.0/L
println("f = $(f)")
println("L = $(L)")
println("g = $(g)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
flush(stdout)
cd(out_path)
mkcd("one_comp_london_L=$(L)_T=$(round(temps[1]; digits=2))-$(round(temps[end]; digits=2))_g=$(round(g; digits=3))_fL=$(round(f*L; digits=1))_test")

# Make ab inito un-correlated phases state
init_syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,f)
cubs = [Cuboid(1, init_syst, split, 1/T_start; u⁺=0.0) for T in temps]



# Time-estimation
################################################################################################

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

#user_input = readline(stdin)
#if user_input == "n"
#    exit()
#end



# Adjusting simulation constants
################################################################################################

println("Adjusting simulation constants.")
@time for cub in cubs; tuneUpdates!(cub); end
println("Controls after adjustment:")
printControls(cubs)



# Cooldown
################################################################################################
println("\nStarted cooldown at: $(Dates.format(now(), "HH:MM"))")

# Do M_th MCS in N_steps steps
M_pr_step = ceil(Int64, M_col/N_steps)
M_col = M_pr_step*N_steps
# Each step is associated with a different set of temperatures given by each row in the matrix
temp_mt = genGeometricTemperatureSteps(T_start, temps, N_steps)
# Energy storage
E_matrix = Array{Float64, 2}(undef, M_col+1, N_T);
for (i, cub) = enumerate(cubs); E_matrix[1,i] = energy(cub); end
E_therm = Array{Float64, 2}(undef, M_th, N_T);
# Vorticity storage
S⁺_col_by_T_aux = [Array{Array{Float64, 2}, 1}(undef, M_pr_step) for i = 1:N_T]
S⁺_col_by_T = [Array{Array{Float64, 2}, 1}(undef, N_steps) for i = 1:N_T]

println("Cooling down from T=$(T_start) using $(N_steps) steps with $(M_pr_step) MCS pr. step.")
flush(stdout)

t_col = @elapsed for step = 1:N_steps
    # First set the temperatures associated with this step to each system.
    β_step = [1/T for T in temp_mt[step, :]]
    for (k, cub) = enumerate(cubs); setTemp!(cub, temp_mt[step, k]); end

    # Do the MCS associated with this step
    for i = 1:M_pr_step
        for (k, cub) = enumerate(cubs)
            δE, _ = mcSweepEnUp!(cub)
            E_matrix[(step-1)*M_pr_step+i+1, k] = E_matrix[(step-1)*M_pr_step+i, k] + δE
            proj_V⁺, proj_V⁻ = xyVortexSnapshot(cub)
            S⁺_col_by_T_aux[k][i], _ = structureFunction(proj_V⁺, proj_V⁻)
        end
    end

    # In the end of a step we correct energy rounding-error
    for (j, cub) = enumerate(cubs); E_matrix[step*M_pr_step+1, j] = energy(cub); end

    # Then we also calculate the average structure function at this step.
    for k = 1:N_T;  S⁺_col_by_T[k][step] = mean(S⁺_col_by_T_aux[k]); end

end

E_matrix = E_matrix./N;

JLD.save("cooldown.jld", "sp", S⁺_col_by_T, "E_matrix", E_matrix, "temp_matrix", temp_mt, "N_steps", N_steps, "M_pr_step", M_pr_step, "M_est", M_est)



# Adjusting simulation constants
################################################################################################

println("Adjusting simulation constants.")
@time for cub in cubs; tuneUpdates!(cub); end
println("Controls after adjustment:")
printControls(cubs)



# Thermalization
################################################################################################
println("Thermalizing at target temperatures for additional $(M_th) MCS")
println("\nStarted thermalization at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

for j = 1:N_T; E_therm[1, j] = energy(cubs[j]); end
t_th = @elapsed for i = 2:M_th
    δE = [mcSweepEnUp!(cub)[1] for cub in cubs]
    E_therm[i, :] .= E_therm[i-1, :] .+ δE
end
AR_est = [estimateAR!(cub)[1] for cub in cubs]

E_therm = E_therm./N

t_th += t_col
M_th += M_col

# Update time-estimation
t_MCS = t_th/(M_th*N_T)

therm_lattices = [getLattice(cub) for cub in cubs]
JLD.save("thermalization.jld", "e_thm", E_therm, "syst", init_syst, "lattices", therm_lattices, "M_th", M_th, "M_est", M_est)
println("Cooldown and thermalization used $(round(t_th/3600; digits=1)) h. Energies and states saved to file.")
println("After thermalization the estimated AR are:")
for est in AR_est; println("AR:\t$(round(est; digits=2))"); end



# Doing measurements
################################################################################################
println("We will now do $(Δt*M) MCS pr temperature  which gives
$(M) measurements and has an ETC of $(round(Δt*M*N_T*t_MCS/3600; digits=2)) h")
println("Started measurements at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

# Setup storage for dual stifness and energy
ρˣ₂_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
ρˣ₃_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
ρʸ₁_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
ρʸ₃_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
ρᶻ₁_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
ρᶻ₂_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
E_by_T = [Array{Float64, 1}(undef, M) for k = 1:N_T]
E_prev = [energy(cub) for cub in cubs]
# Storage for vortices and dual vortices
S⁺_by_T = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_T]
S⁻_by_T = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_T]

# We preform M measurements
@time for m = 1:M
    δE = zeros(N_T)
    for j = 1:Δt
        δE .= δE .+ [mcSweepEnUp!(cub)[1] for cub in cubs]
    end
    all_res = [dualStiffnesses(cub) for cub in cubs]
    # Sort results in terms of temperature and extract ρˣˣₖ₂.
    for (i, res) = enumerate(all_res)
        ρˣ₂_by_T[i][m] = res[1]; ρˣ₃_by_T[i][m] = res[2]
        ρʸ₁_by_T[i][m] = res[3]; ρʸ₃_by_T[i][m] = res[4]
        ρᶻ₁_by_T[i][m] = res[5]; ρᶻ₂_by_T[i][m] = res[6]
        E_by_T[i][m] = E_prev[i] + δE[i]
        E_prev[i] += δE[i]
    end

    # Measure vortices and structure function
    for k = 1:N_T
        proj_V⁺, proj_V⁻ = xyVortexSnapshot(cubs[k])
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        S⁺_by_T[k][m], S⁻_by_T[k][m] = structureFunction(proj_V⁺, proj_V⁻)
    end

end
final_lattices = [getLattice(cub) for cub in cubs]
vortices_by_T = [vortexSnapshot(lattice, getSyst(cubs[1])) for lattice in final_lattices]



# Saving results to file
################################################################################################

JLD.save("energies.jld", "E_by_T", E_by_T)

JLD.save("dual_stiffs.jld", "x2", ρˣ₂_by_T, "x3", ρˣ₃_by_T, "y1", ρʸ₁_by_T, "y3", ρʸ₃_by_T, "z1", ρᶻ₁_by_T, "z2", ρᶻ₂_by_T)
JLD.save("meta.jld", "L1", L₁, "L2", L₂, "L3", L₃, "M", M, "dt", Δt, "temps", temps, "f", f, "kappa", κ₅, "g", g, "nu", ν)
JLD.save("final_states.jld", "lattices", final_lattices)
JLD.save("vorticity.jld", "vortexes", vortices_by_T, "sp", S⁺_by_T, "sm", S⁻_by_T)









