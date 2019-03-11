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
M_th = 2^13  # Number of PT-steps to do for thermalization.
M = 2^12    # Number of measurements
Δt = 16    # Number of PT-steps between each measurement
N_mc = 1   # Number of MCS between each PT-step
# L is assumed to be even.
L = 32     # System length
L₃ = 32
N = L^2*L₃
# Make geometric progression of temperatures between T₁ and Tₘ
T₁ = 0.5
Tₘ = 0.8
N_temp = 3*2^0
R = (Tₘ/T₁)^(1/(N_temp-1))
temps = [T₁*R^(k-1) for k = 1:N_temp]
temps = reverse(temps)
println("Temps.: $(temps)")
κ₅ = 1.0
mkcd("one_comp_london_T=$(round(T₁; digits=2))-$(round(Tₘ; digits=2))_L=$(L)")

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 1.0/L
println("f set to $(f)")
sim = Controls(π-8π/12, 1.0, 0.07)

# Construct k-matrix where the horizontal axis contains kx ∈ [-π, π), while
# the vertical axis contain ky ∈ [-π, π) at the second component
k_matrix = [[2π/L*(x-1-L/2), 2π/L*(L/2-y)] for y=1:L, x=1:L]

# Setup measurement storage
N_T = length(temps)
ρˣ₂_avg_by_T = Array{Float64}(undef, N_T); ρˣ₂_err_by_T = Array{Float64}(undef, N_T);

# Make ab inito un-correlated phases state
init_syst_list = [SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T) for T in temps]
init_ψs = [State(1, syst; u⁺=1.0, u⁻=0.0) for syst in init_syst_list]
init_sim_list = [copy(sim) for syst in init_syst_list];
#ψ_ref = State(2, init_syst_list[1]; u⁺=1.0, u⁻=0.0)

# Making parallel tempering variables
pt = PTRun(init_ψs, init_sim_list, N_mc);


# Time-estimation
################################################################################################

M_est = 2^6
# Do M_est parallel tempering steps
t_meas = @elapsed for i = 1:M_est
    for j = 1:Δt
        PTStep!(pt)
    end

    all_res = dmap(R -> gaugeStiffness([R.ψ]), pt.dr_list)
    En_res = E(pt)
   
    # Measure vortices and structure function
    vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt.dr_list)
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
t_MCS = t_meas/(Δt*pt.N_temp*pt.N_mc)
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS.")

# Estimating ETC
println("Measurements and thermalization will do $((Δt*M + M_th)*pt.N_mc) MCS pr temperature which has an 
ETC of $(round((Δt*M + M_th)*N_T*pt.N_mc*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")

user_input = readline(stdin)
if user_input == "n"
    exit()
end



# Thermalization
################################################################################################
println("Started thermalization at: $(Dates.format(now(), "HH:MM"))")

E_matrix = Array{Float64}(undef, M_th, N_T);
# Do M_th parallel tempering steps
t_th = @elapsed for i = 1:M_th
    PTStep!(pt)
    E_matrix[i, :] = E(pt)
end
E_matrix = E_matrix./N;

# Update time-estimation
t_PT = t_th/M_th
t_MCS = t_PT/(pt.N_temp*pt.N_mc)

int = floor(Int64, M_th/2):M_th
therm_plt = plot(collect(int).+M_est, [E_matrix[int, i] for i = 1:N_T]; 
                 label=reshape(reverse(["T = $(round(T; digits=2))" for T in temps]), (1, N_T)),
                 xaxis="PTS", yaxis="Energy pr. site")
savefig(therm_plt, "thermalization energies.pdf")
println("Thermalization used $(round(t_th/60; digits=1)) m. Energy plots saved to file.")



# Adjusting simulation constants and restarting parallel tempering
################################################################################################

println("Adjusting simulation constants.")
ψs = [R.ψ for R in localize(pt.dr_list)]
adjust_sims = [R.sim for R in localize(pt.dr_list)]
@time adjustSimConstants!(ψs, adjust_sims);
println("Controls after adjustment, choosing number 1 for later simulations:")
printSimControls(adjust_sims)
sim_list = [adjust_sims[N_T] for sim in adjust_sims]
pt = PTRun(ψs, sim_list, N_mc);



# Doing measurements
################################################################################################
println("We will now do $(Δt*M*pt.N_temp*pt.N_mc) MCS in $(Δt*M) PT steps which has an 
ETC of $(round(Δt*M*t_PT/3600; digits=2)) h")
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
    for j = 1:Δt
        PTStep!(pt)
    end
    all_res = dmap(R -> gaugeStiffness([R.ψ]), pt.dr_list)
    En_res = E(pt)
    # Sort results in terms of temperature and extract ρˣˣₖ₂.
    for (i, res) = enumerate(all_res[pt.rep_map])
        ρˣ₂_by_T[i][m] = res[1][1]
        E_by_T[i][m] = En_res[i]
    end

    # Measure vortices and structure function
    vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt.dr_list)
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
vortices_by_T = dmap(R -> vortexSnapshot(R.ψ), pt.dr_list)
ψs = [R.ψ for R in localize(pt.dr_list)]



# Saving results to file
################################################################################################

writedlm("dual_stiffs.data", ρˣ₂_by_T, ':')
writedlm("energies.data", E_by_T, ':')
writedlm("temps.data", temps, ':')

JLD.save("meta.jld", "L", L, "L3", L₃, "M", M, "dt", Δt, "temps", temps, "f", f, "kappa", κ₅, "Psi", ψs)
JLD.save("vorticity.jld", "vortexes", vortices_by_T, "sp", S⁺_by_T, "sm", S⁻_by_T)
#rev_temps = reverse(temps)
#for k = 1:N_T
#    JLD.save("vorti_p_T=$(round(rev_temps[k]; digits=2)).jld", "vor", V⁺_by_T[k])
#end









