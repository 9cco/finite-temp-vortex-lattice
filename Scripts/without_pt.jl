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
g = 0.1    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
M_th = 10^6#2^12  # Number of PT-steps to do for thermalization.
M = 2^10#2^13    # Number of measurements
Δt = 80    # Number of PT-steps between each measurement
N_mc = 1   # Number of MCS between each PT-step
# L is assumed to be even.
L = 12     # System length
L₃ = 12
N = L^2*L₃
# Make geometric progression of temperatures between T₁ and Tₘ
T₁ = 0.2
Tₘ = 0.4
N_temp = 3*2^0
R = (Tₘ/T₁)^(1/(N_temp-1))
temps = [T₁*R^(k-1) for k = 1:N_temp]
temps = reverse(temps)
println("Temps.: $(temps)")
κ₅ = 1.0
mkcd("non_pt_uncorr_one_comp_london_T=$(round(T₁; digits=2))-$(round(Tₘ; digits=2))_L=$(L)_oldKinEn")

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 1.0/L
println("f set to $(f)")
sim = Controls(π-π/12, 1.0, 0.7)

# Construct k-matrix where the horizontal axis contains kx ∈ [-π, π), while
# the vertical axis contain ky ∈ [-π, π) at the second component
k_matrix = [[2π/L*(x-1-L/2), 2π/L*(L/2-y)] for y=1:L, x=1:L]

# Setup measurement storage
N_T = length(temps)
ρˣ₂_avg_by_T = Array{Float64}(undef, N_T); ρˣ₂_err_by_T = Array{Float64}(undef, N_T);

# Make ab inito un-correlated phases state
init_syst_list = [SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T) for T in temps]
init_ψs = [State(6, syst; u⁺=1.0, u⁻=0.0) for syst in init_syst_list]
init_sim_list = [copy(sim) for syst in init_syst_list];
#ψ_ref = State(2, init_syst_list[1]; u⁺=1.0, u⁻=0.0)

# Making parallel tempering variables
#pt = PTRun(init_ψs, init_sim_list, N_mc);
ψ_list = init_ψs
sim_list = init_sim_list



# Time-estimation
################################################################################################

M_est = 2^6
# Do M_est steps of the measurement
t_meas = @elapsed for i = 1:M_est
    nMCS!(ψ_list, sim_list, Δt)
   
    all_res = [gaugeStiffness([ψ]) for ψ in ψ_list]
    En_res = pmap(E, ψ_list)

    # Measure vortices and structure function
    vortices_by_T = [vortexSnapshot(ψ) for ψ in ψ_list]
    #vortices_by_T = pmap(vortexSnapshot, ψ_list)
    for k = 1:N_T
        V⁺, V⁻ = vortices_by_T[k]
        # Take the average over the z-direction.
        proj_V⁺ = avgVort(V⁺); proj_V⁻ = avgVort(V⁻)
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        for ix = 1:L, iy = 1:L
            structureFunction(k_matrix[iy,ix], proj_V⁺, proj_V⁻)
        end
    end
end
t_meas = t_meas/M_est
t_MCS = t_meas/(Δt*N_T)
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS.")

# Estimating ETC
println("Measurements will do $((Δt*M + M_th)*N_T) MCS which has an 
ETC of $(round((Δt*M)*N_T*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")
user_input = readline(stdin)
if user_input == "n"
    exit()
end



# Thermalization
################################################################################################
println("Started thermalization at: $(Dates.format(now(), "HH:MM"))")

# Do max M_th MCS

println("Thermalizing $(N_T) initial states from correlated state")
t_th = @elapsed begin
    thermalized, time, ψ_list, sim_list, E_matrix = @time flatThermalization!(ψ_list, sim_list;
        visible=true, N_SUBS=4, T_AVG=2000, T_QUENCH=100, CUTOFF_MAX=M_th);
end
M_th = size(E_matrix, 2)
E_matrix = E_matrix./N;

# Update time-estimation
t_MCS = t_th/(time*N_T)

int = 1:M_th#floor(Int64, M_th/2):M_th
therm_plt = plot(collect(int).+M_est, [E_matrix[k, int] for k = 1:N_T]; 
                 label=reshape(reverse(["T = $(round(T; digits=2))" for T in temps]), (1, N_T)),
                 xaxis="MCS", yaxis="Energy pr. site")
savefig(therm_plt, "thermalization energies.pdf")
println("Thermalization used $(round(t_th/3600; digits=1)) h. Energy plots saved to file.")



# Adjusting simulation constants and restarting parallel tempering
################################################################################################

println("Adjusting simulation constants.")
@time adjustSimConstants!(ψ_list, sim_list);
println("Controls after adjustment, choosing number 1 for later simulations:")
printSimControls(sim_list)



# Doing measurements
################################################################################################
println("We will now do $(Δt*M) MCS for $(N_T) temperatures which has an 
ETC of $(round(Δt*M*N_T*t_MCS/3600; digits=2)) h")
println("Started measurements at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

# Setup storage for ρˣˣₖ₂
ρˣ₂_by_T = Array{Array{Float64, 1}}(undef, length(temps))
E_by_T = Array{Array{Float64, 1}}(undef, length(temps))
for i = 1:length(temps)
    ρˣ₂_by_T[i] = Array{Float64}(undef, M)
    E_by_T[i] = Array{Float64}(undef, M)
end
# Storage for vortices and dual vortices
S⁺_by_T = [[Array{Float64, 2}(undef, L,L) for m = 1:M] for k = 1:N_T]
S⁻_by_T = [[Array{Float64, 2}(undef, L,L) for m = 1:M] for k = 1:N_T]

# We preform M measurements
@time for m = 1:M
    nMCS!(ψ_list, sim_list, Δt)

    all_res = [gaugeStiffness([ψ]) for ψ in ψ_list]
    En_res = pmap(E, ψ_list)
    # Sort results in terms of temperature and extract ρˣˣₖ₂.
    for (i, res) = enumerate(all_res)
        ρˣ₂_by_T[i][m] = res[1][1]
        E_by_T[i][m] = En_res[i]
    end

    # Measure vortices and structure function
    #vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt.dr_list)
    vortices_by_T = [vortexSnapshot(ψ) for ψ in ψ_list]
    for k = 1:N_T
        V⁺, V⁻ = vortices_by_T[k]
        # Take the average over the z-direction.
        proj_V⁺ = avgVort(V⁺); proj_V⁻ = avgVort(V⁻)
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        for ix = 1:L, iy = 1:L
            S⁺_by_T[k][m][iy,ix], S⁻_by_T[k][m][iy,ix] = structureFunction(k_matrix[iy,ix], proj_V⁺, proj_V⁻)
        end
    end
end
vortices_by_T = pmap(vortexSnapshot, ψ_list)



# Saving results to file
################################################################################################

writedlm("dual_stiffs.data", ρˣ₂_by_T, ':')
writedlm("energies.data", E_by_T, ':')
writedlm("temps.data", temps, ':')

JLD.save("meta.jld", "L", L, "L3", L₃, "M", M, "dt", Δt, "temps", temps, "f", f, "kappa", κ₅, "Psi", ψ_list)
JLD.save("vorticity.jld", "vortexes", vortices_by_T, "sp", S⁺_by_T, "sm", S⁻_by_T)
#rev_temps = reverse(temps)
#for k = 1:N_T
#    JLD.save("vorti_p_T=$(round(rev_temps[k]; digits=2)).jld", "vor", V⁺_by_T[k])
#end









