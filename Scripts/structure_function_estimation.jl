# Current intention:
# Run a simulation with only the 3-4 temperatures and parameters as in [28], but this time with larger system
# size L=64 and also thermalize gradually from a high temperature, but perhaps not as long as last time.


using Distributed
# Script for investigating amplitude dependence of potential
@everywhere using Distributions
@everywhere using Test
@everywhere using StatsBase
@everywhere using BenchmarkTools
@everywhere using LinearAlgebra
@everywhere using LaTeXStrings
using Dates
using Primes
using MCMCDiagnostics
using SharedArrays
using DelimitedFiles
@everywhere src_path = "../Source/"
@everywhere include(src_path*"types.jl")
@everywhere include(src_path*"functions_msc.jl")
@everywhere include(src_path*"functions_neighbors.jl")
@everywhere include(src_path*"functions_types.jl")
@everywhere include(src_path*"functions_energy.jl")
@everywhere include(src_path*"functions_mc.jl")
@everywhere include(src_path*"functions_thermalization.jl")
@everywhere include(src_path*"functions_observables.jl")
@everywhere include(src_path*"functions_symmetries.jl")
@everywhere include(src_path*"jackknife_estimates.jl")
@everywhere include(src_path*"dist_lists.jl")
@everywhere include(src_path*"parallel_tempering.jl")
include(src_path*"functions_plots_and_files.jl")
using Plots
pyplot()
using JLD


# Enter data directory for structure function
fixRC()
# We run a simulation with the parameters
g = √(1/10)    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
M_th = 2^15  # Number of PT-steps to do for thermalization.
M = 2^11    # Number of measurements
Δt = 19    # Number of PT-steps between each measurement
N_mc = 2   # Number of MCS between each PT-step
# L is assumed to be even.
L = 64     # System length
L₃ = 64
N = L^2*L₃
# Make geometric progression of temperatures between T₁ and Tₘ
T₁ = 0.7
Tₘ = 4.0
N_T = 3*2^0
N_steps = 2^10   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_th) must be >= 2
T_start = 2*4.0+0.1   # Start temperature of states when thermalizing. Must be higher than maximum(temps)
R = (Tₘ/T₁)^(1/(N_T-1))
temps = 2*[0.35, 0.4, 1.66]#[T₁*R^(k-1) for k = 1:N_T]#
println("Temps.: $(temps)")
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 2.0/L
println("f set to $(f)")
sim = Controls(π-8π/12, 1.0, 0.1)
mkcd("test_two_comp_london_T=$(round(temps[1]; digits=2))-$(round(temps[end]; digits=2))_L=$(L)_g=$(round(g; digits=3))_fL=$(round(f*L; digits=1))_gradual_high_T_quench")

# Make ab inito un-correlated phases state
init_syst_list = [SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T_start) for T in temps]
init_ψs = [State(6, syst; u⁺=√(0.2), u⁻=√(2)) for syst in init_syst_list]
init_sim_list = [copy(sim) for syst in init_syst_list];
#ψ_ref = State(2, init_syst_list[1]; u⁺=1.0, u⁻=0.0)

# Making parallel tempering variables
pt = PTRun(init_ψs, init_sim_list, N_mc);


# Time-estimation
################################################################################################

M_est = 2^4
# Do M_est parallel tempering steps
t_meas = @elapsed for i = 1:M_est
    for j = 1:Δt
        PTStep!(pt)
    end

    all_res = dmap(R -> gaugeStiffness([R.ψ]), pt)
    En_res = E(pt)
   
    # Measure vortices and structure function
    vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt)
    for k = 1:N_T
        V⁺, V⁻ = vortices_by_T[k]
        # Take the average over the z-direction.
        proj_V⁺ = avgVort(V⁺); proj_V⁻ = avgVort(V⁻)
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        structureFunction(proj_V⁺, proj_V⁻)
    end

end
t_meas = t_meas/M_est
t_MCS = t_meas/(Δt*pt.N_temp*pt.N_mc)
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS.")

# Estimating ETC
println("Measurements and thermalization will do $((Δt*M + M_th)*pt.N_mc) MCS for each of the $(pt.N_temp)
temperatures, which has an ETC of $(round((Δt*M + M_th)*N_T*pt.N_mc*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")
flush(stdout)

user_input = readline(stdin)
if user_input == "n"
    exit()
end



# Thermalization
################################################################################################
println("\nStarted thermalization at: $(Dates.format(now(), "HH:MM"))")

# Do M_th parallel tempering steps in N_steps steps
M_pr_step = ceil(Int64, M_th/(2*N_steps))
M_col = M_pr_step*N_steps
# Each step is associated with a different set of temperatures given by each row in the matrix
temp_mt = genGeometricTemperatureSteps(T_start, temps, N_steps)
E_matrix = Array{Float64}(undef, 2*M_col, N_T);

println("Cooling down from T=$(T_start) using $(N_steps) steps with $(M_pr_step*pt.N_mc) MCS pr. step.")
flush(stdout)

t_col = @elapsed for step = 1:N_steps
    # First set the temperatures associated with this step to each of the replicas in pt
    β_step = [1/T for T in temp_mt[step, :]]
    distributeTemperatures!(pt, β_step)

    # Do the parallel tempering steps associated with this step
    for i = 1:M_pr_step
        PTStep!(pt)
        E_matrix[(step-1)*M_pr_step+i, :] = E(pt)
    end
end

println("Thermalizing at target temperatures for additional $(M_col) PT-steps")
flush(stdout)

t_th = @elapsed for i = 1:M_col
    PTStep!(pt)
    E_matrix[M_col+i, :] = E(pt)
end
t_th += t_col
M_th = 2*M_col

E_matrix = E_matrix./N;

# Update time-estimation
t_PT = t_th/M_th
t_MCS = t_PT/(pt.N_temp*pt.N_mc)

int = 1:M_th#floor(Int64, M_th/2)
therm_plt = plot(collect(int).+M_est, [E_matrix[int, i] for i = 1:N_T]; 
                 label=reshape(["T = $(round(T; digits=2))" for T in temps], (1, N_T)),
                 xaxis="PTS", yaxis="Energy pr. site")
savefig(therm_plt, "thermalization energies.pdf")
JLD.save("therm_energies.jld", "e_mt", E_matrix)
println("Thermalization used $(round(t_th/3600; digits=1)) h. Energies saved to file.")



# Adjusting simulation constants and restarting parallel tempering
################################################################################################

println("Adjusting simulation constants.")
ψs = [R.ψ for R in localize(pt.dr_list)[pt.rep_map]]
adjust_sims = [R.sim for R in localize(pt.dr_list)[pt.rep_map]]
@time adjustSimConstants!(ψs, adjust_sims);
println("Controls after adjustment, choosing number 1 for later simulations:")
printSimControls(adjust_sims)
sim_list = [adjust_sims[N_T] for sim in adjust_sims]
pt = PTRun(ψs, sim_list, N_mc);



# Doing measurements
################################################################################################
println("We will now do $(Δt*M*pt.N_mc) MCS pr temperature in $(Δt*M) PT steps which gives
$(M) measurements and has an ETC of $(round(Δt*M*t_PT/3600; digits=2)) h")
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
# Storage for vortices and dual vortices
S⁺_by_T = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_T]
S⁻_by_T = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_T]

# We preform M measurements
@time for m = 1:M
    for j = 1:Δt
        PTStep!(pt)
    end
    all_res = dmap(R -> gaugeStiffness([R.ψ]), pt)
    En_res = E(pt)
    # Sort results in terms of temperature and extract ρˣˣₖ₂.
    for (i, res) = enumerate(all_res)
        ρˣ₂_by_T[i][m] = res[1][1]; ρˣ₃_by_T[i][m] = res[2][1]
        ρʸ₁_by_T[i][m] = res[3][1]; ρʸ₃_by_T[i][m] = res[4][1]
        ρᶻ₁_by_T[i][m] = res[5][1]; ρᶻ₂_by_T[i][m] = res[6][1]
        E_by_T[i][m] = En_res[i]
    end

    # Measure vortices and structure function
    vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt)
    for k = 1:N_T
        V⁺, V⁻ = vortices_by_T[k]
        # Take the average over the z-direction.
        proj_V⁺ = avgVort(V⁺); proj_V⁻ = avgVort(V⁻)
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        S⁺_by_T[k][m], S⁻_by_T[k][m] = structureFunction(proj_V⁺, proj_V⁻)
    end
end
vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt)
ψs = [R.ψ for R in localize(pt.dr_list)[pt.rep_map]]



# Saving results to file
################################################################################################

writedlm("dual_stiffs.data", ρˣ₂_by_T, ':')
writedlm("energies.data", E_by_T, ':')
writedlm("temps.data", temps, ':')

JLD.save("dual_stiffs.jld", "x2", ρˣ₂_by_T, "x3", ρˣ₃_by_T, "y1", ρʸ₁_by_T, "y3", ρʸ₃_by_T, "z1", ρᶻ₁_by_T, "z2", ρᶻ₂_by_T)
JLD.save("meta.jld", "L", L, "L3", L₃, "M", M, "dt", Δt, "temps", temps, "f", f, "kappa", κ₅, "g", g, "nu", ν, "Psi", ψs)
JLD.save("vorticity.jld", "vortexes", vortices_by_T, "sp", S⁺_by_T, "sm", S⁻_by_T)
#rev_temps = reverse(temps)
#for k = 1:N_T
#    JLD.save("vorti_p_T=$(round(rev_temps[k]; digits=2)).jld", "vor", V⁺_by_T[k])
#end









