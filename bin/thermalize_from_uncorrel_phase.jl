# Script for investigating amplitude dependence of potential

@everywhere using Distributions
@everywhere using Base.Test
@everywhere using StatsBase
@everywhere using BenchmarkTools
using MCMCDiagnostics
@everywhere src_path = "../../../Source/"
@everywhere include(src_path*"types.jl")
@everywhere include(src_path*"functions_msc.jl")
@everywhere include(src_path*"functions_neighbors.jl")
@everywhere include(src_path*"functions_types.jl")
@everywhere include(src_path*"functions_energy.jl")
@everywhere include(src_path*"functions_mc.jl")
@everywhere include(src_path*"functions_thermalization.jl")
@everywhere include(src_path*"functions_observables.jl")
include(src_path*"functions_plots_and_files.jl")
using Plots
gr()


THERM_FRAC = 1/10
DT_MAX = 10000
MEASURE_FILE = "measured.statelist"

@everywhere const two_pi = 2π

# We run a simulation with the parameters
g = 0.3    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
L = 24     # System length
L₃ = 24
T_list = [T for T = 0.01:0.1:2]
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 0.0/L
println("f set to $(f)")
sim = Controls(π-π/12, 1.0, 4.0)

M = 50
# Setup measurement storage
N_T = length(T_list)
u⁺_avg_by_T = Array{Float64}(N_T); u⁻_avg_by_T = Array{Float64}(N_T)
u⁺_err_by_T = Array{Float64}(N_T); u⁻_err_by_T = Array{Float64}(N_T)

# Make ab inito un-correlated phases state
nw = nprocs()-1 # The number of low temperature states needed
init_syst_list = Array{SystConstants, 1}(nw+1)
init_syst_list[1:nw] = [SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T_list[1]) for i = 1:nw] # Make nw low temp-states
init_syst_list[nw+1] = SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T_list[end])              # Make final high temp state.
init_ψ_list = [State(5, syst; u⁺=1.0, u⁻=0.0) for syst in init_syst_list]
init_sim_list = [copy(sim) for syst in init_syst_list];


# Run these states until the energy-curve is flat.

println("Thermalizing $(nw+1) initial states from correlated state")
thermalized, time, init_ψ_list, init_sim_list, E_matrix = @time flatThermalization!(init_ψ_list, init_sim_list; visible=true);



# Plot last energies
plt = plot(1:size(E_matrix,2), [E_matrix[s,:] for (s, syst) in enumerate(init_syst_list)])
savefig(plt, "initial_states_thermalization_last_average.pdf")



# Then preform thermalization and measurement for each temperature individually in separate folders
@time for (i, T) in enumerate(T_list)
    println("\nEntering simulation of T = $(T)")
    dir_name = "T_$(T)"
    mkcd(dir_name)

    # Start with fresh states
    ψ_ref = copy(init_ψ_list[end])
    ψ_w = [copy(ψ) for ψ in init_ψ_list[1:nw]]
    # Update to new temperature
    syst = SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T)
    ψ_ref.consts = syst
    for ψ in ψ_w
        ψ.consts = syst
    end

    # Thermalize states
    println("Thermalizing from high and low temperature")
    thermalized, t₀, ψ_ref, E_ref, sim_ref, ψ_w, E_w, sim_w = thermalizeLite!(ψ_ref, ψ_w, copy(sim))

    # Preform measurements by saving states to file.
    ψ_list = [ψ_ref, ψ_w...]
    Δt = min(DT_MAX, ceil(Int64, t₀*THERM_FRAC))
    measurementSeries!(ψ_list, sim_ref, M, Δt; filename=MEASURE_FILE)

    # Load states to memory.
    ψ_measured = loadStates(MEASURE_FILE)

    # Use states to produce lists of order-parameters
    #η⁺_list, η⁻_list = measureOrdParam(ψ_measured)
    # Get averages and errors
    #η⁺_avg, η⁺_err = avgErr(η⁺_list); η⁻_avg, η⁻_err = avgErr(η⁻_list)
    # Get correlation times
    #η⁺_τ = M/effective_sample_size(η⁺_list); η⁻_τ = M/effective_sample_size(η⁻_list)

    # Use states to produce lists of order-parameter amplitudes
    u⁺_list, u⁻_list = measureMeanAmplitudes(ψ_measured)
    # Get averages and errors
    u⁺_avg, u⁺_err = avgErr(u⁺_list); u⁻_avg, u⁻_err = avgErr(u⁻_list)
    # Get correlation times
    u⁺_τ = M/effective_sample_size(u⁺_list); u⁻_τ = M/effective_sample_size(u⁻_list)
    u⁺_avg_r, u⁺_err_r = scientificRounding(u⁺_avg, u⁺_err); u⁻_avg_r, u⁻_err_r = scientificRounding(u⁻_avg, u⁻_err)
    println("⟨|η⁺|⟩ = $(u⁺_avg_r) ± $(u⁺_err_r)\t\tτ = $(u⁺_τ)")
    println("⟨|η⁻|⟩ = $(u⁻_avg_r) ± $(u⁻_err_r)\t\tτ = $(u⁻_τ)")

    # Save result to array
    u⁺_avg_by_T[i] = u⁺_avg; u⁻_avg_by_T[i] = u⁻_avg;
    u⁺_err_by_T[i] = u⁺_err; u⁻_err_by_T[i] = u⁻_err;

    cd("../")
end


# Plot results

plt = scatter(T_list, u⁺_avg_by_T, yerror=u⁺_err_by_T)
savefig(plt, "+amplitude_by_T.pdf")

plt = scatter(T_list, u⁻_avg_by_T, yerror=u⁻_err_by_T)
savefig(plt, "-amplitude_by_T.pdf")
