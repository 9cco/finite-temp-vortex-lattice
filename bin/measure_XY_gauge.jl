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
pyplot()



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

ρˣ₂_avg_by_T = Array{Float64}(N_T); ρˣ₃_avg_by_T = Array{Float64}(N_T); ρʸ₁_avg_by_T = Array{Float64}(N_T)
ρˣ₂_err_by_T = Array{Float64}(N_T); ρˣ₃_err_by_T = Array{Float64}(N_T); ρʸ₁_err_by_T = Array{Float64}(N_T)
ρʸ₃_avg_by_T = Array{Float64}(N_T); ρᶻ₁_avg_by_T = Array{Float64}(N_T); ρᶻ₂_avg_by_T = Array{Float64}(N_T)
ρʸ₃_err_by_T = Array{Float64}(N_T); ρᶻ₁_err_by_T = Array{Float64}(N_T); ρᶻ₂_err_by_T = Array{Float64}(N_T)




# Then preform thermalization and measurement for each temperature individually in separate folders
@time for (i, T) in enumerate(T_list)
    println("\n\n*****************************************************************************\nEntering simulation of T = $(T)")
    flush(STDOUT)
    dir_name = "T_$(T)"
    mkcd(dir_name)

    # Load states to memory.
    ψ_measured = loadStates(MEASURE_FILE)
    M = length(ψ_measured)

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
    ρˣ₂_list, ρˣ₃_list, ρʸ₁_list, ρʸ₃_list, ρᶻ₁_list, ρᶻ₂_list = gaugeStiffness(ψ_measured)
    # Get averages and errors
    ρˣ₂_avg, ρˣ₂_err = avgErr(ρˣ₂_list); ρˣ₃_avg, ρˣ₃_err = avgErr(ρˣ₃_list); ρʸ₁_avg, ρʸ₁_err = avgErr(ρʸ₁_list)
    ρʸ₃_avg, ρʸ₃_err = avgErr(ρʸ₃_list); ρᶻ₁_avg, ρᶻ₁_err = avgErr(ρᶻ₁_list); ρᶻ₂_avg, ρᶻ₂_err = avgErr(ρᶻ₂_list)
    # Get correlation times and rounding for print
    ρˣ₂_τ = M/effective_sample_size(ρˣ₂_list)
    ρˣ₃_τ = M/effective_sample_size(ρˣ₃_list)
    ρʸ₁_τ = M/effective_sample_size(ρʸ₁_list)
    ρʸ₃_τ = M/effective_sample_size(ρʸ₃_list)
    ρᶻ₁_τ = M/effective_sample_size(ρᶻ₁_list)
    ρᶻ₂_τ = M/effective_sample_size(ρᶻ₂_list)
    ρˣ₂_avg_r, ρˣ₂_err_r = scientificRounding(ρˣ₂_avg, ρˣ₂_err);
    ρˣ₃_avg_r, ρˣ₃_err_r = scientificRounding(ρˣ₃_avg, ρˣ₃_err);
    ρʸ₁_avg_r, ρʸ₁_err_r = scientificRounding(ρʸ₁_avg, ρʸ₁_err);
    ρʸ₃_avg_r, ρʸ₃_err_r = scientificRounding(ρʸ₃_avg, ρʸ₃_err);
    ρᶻ₁_avg_r, ρᶻ₁_err_r = scientificRounding(ρᶻ₁_avg, ρᶻ₁_err);
    ρᶻ₂_avg_r, ρᶻ₂_err_r = scientificRounding(ρᶻ₂_avg, ρᶻ₂_err);
    println("ρˣ(k₂) = $(ρˣ₂_avg_r) ± $(ρˣ₂_err_r)\t\tτ = $(ρˣ₂_τ)")
    println("ρˣ(k₃) = $(ρˣ₃_avg_r) ± $(ρˣ₃_err_r)\t\tτ = $(ρˣ₃_τ)")
    println("ρʸ(k₁) = $(ρʸ₁_avg_r) ± $(ρʸ₁_err_r)\t\tτ = $(ρʸ₁_τ)")
    println("ρʸ(k₃) = $(ρʸ₃_avg_r) ± $(ρʸ₃_err_r)\t\tτ = $(ρʸ₃_τ)")
    println("ρᶻ(k₁) = $(ρᶻ₁_avg_r) ± $(ρᶻ₁_err_r)\t\tτ = $(ρᶻ₁_τ)")
    println("ρᶻ(k₂) = $(ρᶻ₂_avg_r) ± $(ρᶻ₂_err_r)\t\tτ = $(ρᶻ₂_τ)")

    # Save results to arrays
    ρˣ₂_avg_by_T[i] = ρˣ₂_avg; ρˣ₃_avg_by_T[i] = ρˣ₃_avg; ρʸ₁_avg_by_T[i] = ρʸ₁_avg
    ρˣ₂_err_by_T[i] = ρˣ₂_err; ρˣ₃_err_by_T[i] = ρˣ₃_err; ρʸ₁_err_by_T[i] = ρʸ₁_err
    ρʸ₃_avg_by_T[i] = ρʸ₃_avg; ρᶻ₁_avg_by_T[i] = ρᶻ₁_avg; ρᶻ₂_avg_by_T[i] = ρᶻ₂_avg
    ρʸ₃_err_by_T[i] = ρʸ₃_err; ρᶻ₁_err_by_T[i] = ρᶻ₁_err; ρᶻ₂_err_by_T[i] = ρᶻ₂_err

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
    
    # Delete states from memory
    ψ_measured = 0
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

plt = scatter(T_list, L.*ρˣ₂_avg_by_T, yerror=L.*ρˣ₂_err_by_T; xlabel="T", ylabel="L×⟨ρˣ₂⟩", title="|Σᵣ(∇×A)ˣeⁱᵏʳ)|²/((2π)²L³) for k = (0,2π/L,0)");
savefig(plt, "x2_gauge_stiff_by_T.pdf")
plt = scatter(T_list, L.*ρˣ₃_avg_by_T, yerror=L.*ρˣ₃_err_by_T; xlabel="T", ylabel="L×⟨ρˣ₃⟩", title="|Σᵣ(∇×A)ˣeⁱᵏʳ)|²/((2π)²L³) for k = (0,0,2π/L)");
savefig(plt, "x3_gauge_stiff_by_T.pdf")
plt = scatter(T_list, L.*ρʸ₁_avg_by_T, yerror=L.*ρʸ₁_err_by_T; xlabel="T", ylabel="L×⟨ρʸ₁⟩", title="|Σᵣ(∇×A)ʸeⁱᵏʳ)|²/((2π)²L³) for k = (2π/L,0,0)");
savefig(plt, "y1_gauge_stiff_by_T.pdf")
plt = scatter(T_list, L.*ρʸ₃_avg_by_T, yerror=L.*ρʸ₃_err_by_T; xlabel="T", ylabel="L×⟨ρʸ₃⟩", title="|Σᵣ(∇×A)ʸeⁱᵏʳ)|²/((2π)²L³) for k = (0,0,2π/L)");
savefig(plt, "y3_gauge_stiff_by_T.pdf")
plt = scatter(T_list, L.*ρᶻ₁_avg_by_T, yerror=L.*ρᶻ₁_err_by_T; xlabel="T", ylabel="L×⟨ρᶻ₁⟩", title="|Σᵣ(∇×A)ᶻeⁱᵏʳ)|²/((2π)²L³) for k = (2π/L,0,0)");
savefig(plt, "z1_gauge_stiff_by_T.pdf")
plt = scatter(T_list, L.*ρᶻ₂_avg_by_T, yerror=L.*ρᶻ₂_err_by_T; xlabel="T", ylabel="L×⟨ρᶻ₂⟩", title="|Σᵣ(∇×A)ᶻeⁱᵏʳ)|²/((2π)²L³) for k = (0,2π/L,0)");
savefig(plt, "z2_gauge_stiff_by_T.pdf")


# Save data
mkcd("Data_Arrays")
#writedlm("Temp_list.data", T_list, ":")
#writedlm("eta+_avg.data", η⁺_avg_by_T, ":")
#writedlm("eta+_err.data", η⁺_err_by_T, ":")
#writedlm("eta-_avg.data", η⁻_avg_by_T, ":")
#writedlm("eta-_err.data", η⁻_err_by_T, ":")
#writedlm("hel_mod_avg.data", Υ_avg_by_T, ":")
#writedlm("hel_mod_err.data", Υ_err_by_T, ":")
writedlm("x2_gs_avg.data", ρˣ₂_avg_by_T, ":")
writedlm("x2_gs_err.data", ρˣ₂_err_by_T, ":")
writedlm("x3_gs_avg.data", ρˣ₃_avg_by_T, ":")
writedlm("x3_gs_err.data", ρˣ₃_err_by_T, ":")
writedlm("y1_gs_avg.data", ρʸ₁_avg_by_T, ":")
writedlm("y1_gs_err.data", ρʸ₁_err_by_T, ":")
cd("..")

