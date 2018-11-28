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
pyplot()


MEASURE_FILE = "measured.statelist"

@everywhere const two_pi = 2π

# We run a simulation with the parameters
g = 0.3    # Gauge coupling
ν = 0.3    # Anisotropy

# Other parameters
L = 12     # System length
L₃ = 12
T_list = [T for T = 0.01:0.1:2]
κ₅ = 1.0

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = 0.0/L

# Setup measurement storage
N_T = length(T_list)

η⁺_avg_by_T = Array{Float64}(N_T); η⁻_avg_by_T = Array{Float64}(N_T)
η⁺_err_by_T = Array{Float64}(N_T); η⁻_err_by_T = Array{Float64}(N_T)


# For each temperature calculate plot of the energy with amplitude and where the state is on this plot
# as well as a histogram of amplitudes. Also calculate the order parameter
@time for (i, T) in enumerate(T_list)
    println("\nEntering simulation of T = $(T)")
    dir_name = "T_$(T)"
    cd(dir_name)

    syst = SystConstants(L, L₃, 1/g^2, ν, κ₅, f, 1/T)

    # Load states to memory.
    ψ_measured = loadStates(MEASURE_FILE)
    M = length(ψ_measured)

    # Use states to produce lists of order-parameters
    η⁺_list, η⁻_list = measureOrdParam(ψ_measured)
    # Get averages and errors
    η⁺_avg, η⁺_err = avgErr(η⁺_list); η⁻_avg, η⁻_err = avgErr(η⁻_list)
    # Get correlation times
    η⁺_τ = M/effective_sample_size(η⁺_list); η⁻_τ = M/effective_sample_size(η⁻_list)
    η⁺_avg_r, η⁺_err_r = scientificRounding(η⁺_avg, η⁺_err); η⁻_avg_r, η⁻_err_r = scientificRounding(η⁻_avg, η⁻_err)
    println("|⟨η⁺⟩| = $(η⁺_avg_r) ± $(η⁺_err_r)\t\tτ = $(η⁺_τ)")
    println("|⟨η⁻⟩| = $(η⁻_avg_r) ± $(η⁻_err_r)\t\tτ = $(η⁻_τ)")

    # Save result to array
    η⁺_avg_by_T[i] = η⁺_avg; η⁻_avg_by_T[i] = η⁻_avg;
    η⁺_err_by_T[i] = η⁺_err; η⁻_err_by_T[i] = η⁻_err;

    # Find a large number of amplitudes
    n = min(3, length(ψ_measured)) # Number of states to use
    N = length(ψ_measured[1].lattice)
    u⁺_list = Array{Float64}(n*N); u⁻_list = Array{Float64}(n*N);
    i = 1 # Next available index in u-lists.
    for s = 1:n
        for ϕ in ψ_measured[s].lattice
            u⁺_list[i] = ϕ.u⁺; u⁻_list[i] = ϕ.u⁻
            i += 1
        end
    end

    plt = histogram(u⁺_list; title="Amplitudes of + component", yaxis="#", xaxis="|u⁺|", label="T=$(T)")
    savefig(plt, "+comp_histogram.pdf")
    plt = histogram(u⁻_list; title="Amplitudes of - component", yaxis="#", xaxis="|u⁻|", label="T=$(T)")
    savefig(plt, "-comp_histogram.pdf")

    # Now we make a plot of the energy-graph and where the average amplitude is
    u⁺_list, _ = measureMeanAmplitudes(ψ_measured)
    u⁺_avg, u⁺_err = avgErr(u⁺_list)

    u_int = [x for x = 0:0.05:1.5]
    E_int = Array{Float64}(length(u_int))
    for (i,u) in enumerate(u_int)
        # Create correlated state with all amplitudes uniform
        ψ = State(1, syst; u⁺=u)
        E_int[i] = E(ψ)/length(ψ.lattice)
    end
    plt = plot(u_int, E_int; title="Energy of numerical model", xlabel="u⁺", ylabel="E", ylims=(min(minimum(E_int),0)-0.03,maximum(E_int)), label="Potential at T=$(T)")

    ψ = State(1, syst; u⁺=u⁺_avg)
    E_u = E(ψ)/length(ψ.lattice)
    scatter!(plt, [u⁺_avg], [E_u]; xerror=[u⁺_err], ms=7, label="Average energy")
    savefig(plt, "energy_amplitude.pdf")

    cd("../")
end


# Plot results

plt = scatter(T_list, η⁺_avg_by_T, yerror=η⁺_err_by_T; ylims=(0, maximum(η⁺_avg_by_T)+maximum(η⁺_err_by_T)), ylabel="|⟨η⁺⟩|", xlabel="T", label="Only potential")
savefig(plt, "eta+_by_T.pdf")
