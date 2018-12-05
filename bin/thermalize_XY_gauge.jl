# Script for investigating amplitude dependence of potential

@everywhere using Distributions
@everywhere using Base.Test
@everywhere using StatsBase
@everywhere using BenchmarkTools
using MCMCDiagnostics
@everywhere src_path = "../../Source/"
@everywhere include(src_path*"types.jl")
@everywhere include(src_path*"functions_msc.jl")
@everywhere include(src_path*"functions_neighbors.jl")
@everywhere include(src_path*"functions_types.jl")
@everywhere include(src_path*"functions_energy.jl")
@everywhere include(src_path*"functions_mc.jl")
@everywhere include(src_path*"functions_thermalization.jl")
@everywhere include(src_path*"functions_observables.jl")
@everywhere include(src_path*"functions_symmetries.jl")
include(src_path*"functions_plots_and_files.jl")
using Plots
gr()



THERM_FRAC = 1/10
DT_MAX = 10000
MEASURE_FILE = "measured.statelist"

@everywhere const two_pi = 2π

# We run a simulation with the parameters
g = 1.0    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
L = 20     # System length
L₃ = 20
T_list = [T for T = 0.1:0.05:1.8]
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 0.0/L
println("f set to $(f)")
sim = Controls(π-π/12, 1.0, 4.0)

M = 300
# Setup measurement storage
N_T = length(T_list)
#u⁺_avg_by_T = Array{Float64}(N_T); u⁻_avg_by_T = Array{Float64}(N_T)
#u⁺_err_by_T = Array{Float64}(N_T); u⁻_err_by_T = Array{Float64}(N_T)

#η⁺_avg_by_T = Array{Float64}(N_T); η⁻_avg_by_T = Array{Float64}(N_T)
#η⁺_err_by_T = Array{Float64}(N_T); η⁻_err_by_T = Array{Float64}(N_T)

#Υ_avg_by_T = Array{Float64}(N_T)
#Υ_err_by_T = Array{Float64}(N_T);

#ρˣ₂_avg_by_T = Array{Float64}(N_T); ρˣ₃_avg_by_T = Array{Float64}(N_T); ρʸ₁_avg_by_T = Array{Float64}(N_T)
#ρˣ₂_err_by_T = Array{Float64}(N_T); ρˣ₃_err_by_T = Array{Float64}(N_T); ρʸ₁_err_by_T = Array{Float64}(N_T)

# Make ab inito un-correlated phases state
nw = max(1,nprocs()-1) # The number of low temperature states needed
init_syst_list = Array{SystConstants, 1}(nw+1)
init_syst_list[1:nw] = [SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/0.1) for i = 1:nw] # Make nw low temp-states
init_syst_list[nw+1] = SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/0.5)              # Make final high temp state.
init_ψ_list = [State(1, syst; u⁺=1.0, u⁻=0.0) for syst in init_syst_list]
init_sim_list = [copy(sim) for syst in init_syst_list];


# Run these states until the energy-curve is flat.

println("Thermalizing $(nw+1) initial states from correlated state")
thermalized, time, init_ψ_list, init_sim_list, E_matrix = @time flatThermalization!(init_ψ_list, init_sim_list;
    visible=true, N_SUBS=10, T_AVG=5000, T_QUENCH=1000);

N = length(init_ψ_list[1].lattice)
# Plot last energies
plt = plot(1:size(E_matrix,2), [E_matrix[s,:]./N for s = 1:size(E_matrix,1)]; xlabel="MCS", ylabel="Total En", title="Last thermalization interval")
savefig(plt, "initial_states_thermalization_last_average.pdf")

println("Initial energies:")
println([E(init_ψ_list[i])/N for i = 1:length(init_ψ_list)])


# Then preform thermalization and measurement for each temperature individually in separate folders
@time for (i, T) in enumerate(T_list)
    println("\n\n*****************************************************************************\nEntering simulation of T = $(T)")
    flush(STDOUT)
    dir_name = "T_$(T)"
    mkcd(dir_name)

    # Start with fresh states
    ψ_ref = copy(init_ψ_list[end])
    ψ_w = [copy(ψ) for ψ in init_ψ_list[1:nw]]
    println("Thermalizing from high (T=$(1/ψ_ref.consts.β)) and low temperature (T=$(1/ψ_w[1].consts.β))")
    # Update to new temperature
    syst = SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T)
    ψ_ref.consts = syst
    for ψ in ψ_w
        ψ.consts = syst
    end

    # Thermalize states
    @time thermalized, t₀, ψ_ref, E_ref, sim_ref, ψ_w, E_w, sim_w = thermalizeLite!(ψ_ref, ψ_w, copy(sim); visible=true,
        STABILITY_CUTOFF=30000)
    # If we didn't manage to thermalize, try to go back to larger simulation constants
    #if !thermalized
    #    ψ_list = [ψ_ref, ψ_w...]
    #    sim_list = [copy(sim) for ψ in ψ_list]
    #    # Run for 1/4 of the previous thermalization time with larger constants
    #    ψ_list, E_matrix = nMCSEnergy(ψ_list, sim_list, floor(Int64, t₀/4)+1, [E(ψ) for ψ in ψ_list])

    println("Thermalized energies:")
    println([E(ψ_ref), [E(ψ) for ψ in ψ_w]...])

    # Preform measurements by saving states to file.
    ψ_list = [ψ_ref, ψ_w...]
    Δt = min(DT_MAX, ceil(Int64, t₀*THERM_FRAC))
    println("Δt = $(Δt), which means that we will do in total $(M*Δt) MCS")
    @time measurementSeries!(ψ_list, sim_ref, M, Δt; filename=MEASURE_FILE)

    # Load states to memory.
    #ψ_measured = loadStates(MEASURE_FILE)
    #M = length(ψ_measured)

    # Use states to produce lists of order-parameters
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #η⁺_list, η⁻_list = measureOrdParam(ψ_measured)
    # Get averages and errors
    #η⁺_avg, η⁺_err = avgErr(η⁺_list); η⁻_avg, η⁻_err = avgErr(η⁻_list)
    # Get correlation times
    #η⁺_τ = M/effective_sample_size(η⁺_list); η⁻_τ = M/effective_sample_size(η⁻_list)
    #η⁺_avg_r, η⁺_err_r = scientificRounding(η⁺_avg, η⁺_err); η⁻_avg_r, η⁻_err_r = scientificRounding(η⁻_avg, η⁻_err)
    #println("|⟨η⁺⟩| = $(η⁺_avg_r) ± $(η⁺_err_r)\t\tτ = $(η⁺_τ)")
    #println("|⟨η⁻⟩| = $(η⁻_avg_r) ± $(η⁻_err_r)\t\tτ = $(η⁻_τ)")
    
    # Save result to array
    #η⁺_avg_by_T[i] = η⁺_avg; η⁻_avg_by_T[i] = η⁻_avg;
    #η⁺_err_by_T[i] = η⁺_err; η⁻_err_by_T[i] = η⁻_err;
    
    # Use states to produce list of helicity moduli
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #Υ_list = measureHelicityModulus(ψ_measured)
    # Get averages and errors
    #Υ_avg, Υ_err = avgErr(Υ_list)
    # Get correlation times
    #Υ_τ = M/effective_sample_size(Υ_list)
    #Υ_avg_r, Υ_err_r = scientificRounding(Υ_avg, Υ_err)
    #println("⟨Υ⟩ = $(Υ_avg_r) ± $(Υ_err_r)\t\tτ = $(Υ_τ)")
    
    # Save result to array
    #Υ_avg_by_T[i] = Υ_avg; Υ_err_by_T[i] = Υ_err

    # Use states to produce lists of gauge-stiffnesses
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #ρˣ₂_list, ρˣ₃_list, ρʸ₁_list, ρʸ₃_list, ρᶻ₁_list, ρᶻ₂_list = gaugeStiffness(ψ_measured)
    # Get averages and errors
    #ρˣ₂_avg, ρˣ₂_err = avgErr(ρˣ₂_list); ρˣ₃_avg, ρˣ₃_err = avgErr(ρˣ₃_list); ρʸ₁_avg, ρʸ₁_err = avgErr(ρʸ₁_list)
    # Get correlation times and rounding for print
    #ρˣ₂_τ = M/effective_sample_size(ρˣ₂_list)
    #ρˣ₂_avg_r, ρˣ₂_err_r = scientificRounding(ρˣ₂_avg, ρˣ₂_err);
    #println("ρˣ(k₂) = $(ρˣ₂_avg_r) ± $(ρˣ₂_err_r)\t\tτ = $(ρˣ₂_τ)")

    # Save results to arrays
    #ρˣ₂_avg_by_T[i] = ρˣ₂_avg; ρˣ₃_avg_by_T[i] = ρˣ₃_avg; ρʸ₁_avg_by_T[i] = ρʸ₁_avg
    #ρˣ₂_err_by_T[i] = ρˣ₂_err; ρˣ₃_err_by_T[i] = ρˣ₃_err; ρʸ₁_err_by_T[i] = ρʸ₁_err

    # Use states to produce lists of order-parameter amplitudes
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #u⁺_list, u⁻_list = measureMeanAmplitudes(ψ_measured)
    # Get averages and errors
    #u⁺_avg, u⁺_err = avgErr(u⁺_list); u⁻_avg, u⁻_err = avgErr(u⁻_list)
    # Get correlation times
    #u⁺_τ = M/effective_sample_size(u⁺_list); u⁻_τ = M/effective_sample_size(u⁻_list)
    #u⁺_avg_r, u⁺_err_r = scientificRounding(u⁺_avg, u⁺_err); u⁻_avg_r, u⁻_err_r = scientificRounding(u⁻_avg, u⁻_err)
    #println("⟨|η⁺|⟩ = $(u⁺_avg_r) ± $(u⁺_err_r)\t\tτ = $(u⁺_τ)")
    #println("⟨|η⁻|⟩ = $(u⁻_avg_r) ± $(u⁻_err_r)\t\tτ = $(u⁻_τ)")

    # Save result to array
    #u⁺_avg_by_T[i] = u⁺_avg; u⁻_avg_by_T[i] = u⁻_avg;
    #u⁺_err_by_T[i] = u⁺_err; u⁻_err_by_T[i] = u⁻_err;

    @everywhere gc()

    cd("../")
end


# Plot results

#plt = scatter(T_list, u⁺_avg_by_T, yerror=u⁺_err_by_T; xlabel="T", ylabel="⟨|η⁺|⟩", title="Amplitude T dependence")
#savefig(plt, "+amplitude_by_T.pdf")
#plt = scatter(T_list, u⁻_avg_by_T, yerror=u⁻_err_by_T; xlabel="T", ylabel="⟨|η⁺|⟩", title="Amplitude T dependence")
#savefig(plt, "-amplitude_by_T.pdf")

#plt = scatter(T_list, η⁺_avg_by_T, yerror=η⁺_err_by_T; xlabel="T", ylabel="|⟨η⁺⟩|", title="OP T dependence")
#savefig(plt, "OP+_by_T.pdf")
#plt = scatter(T_list, η⁻_avg_by_T, yerror=η⁻_err_by_T; xlabel="T", ylabel="|⟨η⁻⟩|", title="OP T dependence")
#savefig(plt, "OP-_by_T.pdf")

#plt = scatter(T_list, Υ_avg_by_T, yerror=Υ_err_by_T; xlabel="T", ylabel="Υ", title="Helicity modulus T dependence")
#savefig(plt, "Hel_Mod_by_T.pdf")

#plt = scatter(T_list, ρˣ₂_avg_by_T, yerror=ρˣ₂_err_by_T; xlabel="T", ylabel="⟨ρˣ₂⟩", title="|Σᵣ(∇×A)ˣeⁱᵏʳ)| for k = (0,2π/L,0)")
#savefig(plt, "x2_gauge_stiff_by_T.pdf")
#plt = scatter(T_list, ρˣ₃_avg_by_T, yerror=ρˣ₃_err_by_T; xlabel="T", ylabel="⟨ρˣ₃⟩", title="|Σᵣ(∇×A)ˣeⁱᵏʳ)| for k = (0,0,2π/L)")
#savefig(plt, "x3_gauge_stiff_by_T.pdf")
#plt = scatter(T_list, ρʸ₁_avg_by_T, yerror=ρʸ₁_err_by_T; xlabel="T", ylabel="⟨ρʸ₁⟩", title="|Σᵣ(∇×A)ʸeⁱᵏʳ)| for k = (2π/L,0,0)")
#savefig(plt, "x3_gauge_stiff_by_T.pdf")


# Save data
#mkcd("Data_Arrays")
#writedlm("Temp_list.data", T_list, ":")
#writedlm("eta+_avg.data", η⁺_avg_by_T, ":")
#writedlm("eta+_err.data", η⁺_err_by_T, ":")
#writedlm("eta-_avg.data", η⁻_avg_by_T, ":")
#writedlm("eta-_err.data", η⁻_err_by_T, ":")
#writedlm("hel_mod_avg.data", Υ_avg_by_T, ":")
#writedlm("hel_mod_err.data", Υ_err_by_T, ":")
#writedlm("x2_gs_avg.data", ρˣ₂_avg_by_T, ":")
#writedlm("x2_gs_err.data", ρˣ₂_err_by_T, ":")
#writedlm("x3_gs_avg.data", ρˣ₃_avg_by_T, ":")
#writedlm("x3_gs_err.data", ρˣ₃_err_by_T, ":")
#writedlm("y1_gs_avg.data", ρʸ₁_avg_by_T, ":")
#writedlm("y1_gs_err.data", ρʸ₁_err_by_T, ":")
#cd("..")

