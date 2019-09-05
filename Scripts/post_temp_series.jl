using Distributed
# Script for investigating amplitude dependence of potential
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using LinearAlgebra
@everywhere using LaTeXStrings
using BenchmarkTools
using Test
using Dates
using Primes
using MCMCDiagnostics
using SharedArrays
using DelimitedFiles

@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()

src_path = "/home/subr00t/mc/Source/Grid/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere src_path = "/home/subr00t/mc/Source/Grid/"
@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")
include(src_path*"plot_functions.jl")

@everywhere src_path = "/home/subr00t/mc/Source/"
@everywhere include(src_path*"jackknife_estimates.jl")

using Plots
pyplot()

using JLD

plot_folder = "Temp_Series"

folders = filter(s -> isdir(s) && s != plot_folder, readdir())
num_blocks = 2^7
N_F = length(folders)

# Plot storage
Cv_avg_by_T = Array{Float64, 1}(undef, N_F)
Cv_err_by_T = Array{Float64, 1}(undef, N_F)
hex_point_T = Array{Float64, 1}(undef, N_F)
hex_err_T = Array{Float64, 1}(undef, N_F)
quad_point_T = Array{Float64, 1}(undef, N_F)
quad_err_T = Array{Float64, 1}(undef, N_F)
stripe_point_T = Array{Float64, 1}(undef, N_F)
stripe_err_T = Array{Float64, 1}(undef, N_F)
control_T = Array{Float64, 1}(undef, N_F)
control_err_T = Array{Float64, 1}(undef, N_F)
T_list = Array{Float64, 1}(undef, N_F)
u⁺_by_T = Array{Float64, 1}(undef, N_F)
u⁻_by_T = Array{Float64, 1}(undef, N_F)
E_avg_by_T = Array{Float64, 1}(undef, N_F)
E_err_by_T = Array{Float64, 1}(undef, N_F)

# Enter first folder and collect metadata
cd(folders[1])

meta_di = JLD.load("meta.jld")
f = meta_di["f"]; L₃ = meta_di["L3"]; M = meta_di["M"]; κ₅ = meta_di["kappa"];
L₁ = meta_di["L1"]; L₂ = meta_di["L2"]; g = meta_di["g"]; ν = meta_di["nu"]
N = L₁*L₂*L₃;

cd("../")

for (i_f, folder) = enumerate(folders)
    cd(folder)

    # Calculating specific heat

    ens_di = JLD.load("energies.jld")
    E_by_T = ens_di["E_by_T"]

    meta_di = JLD.load("meta.jld")
    f_new = meta_di["f"]; L₃_new = meta_di["L3"]; M_new = meta_di["M"]; temps = meta_di["temps"]; κ₅_new = meta_di["kappa"];
    Δt = meta_di["dt"]; L₁_new = meta_di["L1"]; L₂_new = meta_di["L2"]; g_new = meta_di["g"]; ν_new = meta_di["nu"]
    M_amp = meta_di["M_amp"];
    println("Calculating Cv and structure at T = $(temps[1])")

    N = L₁*L₂*L₃;
    N_T = length(temps)

    if N_T != 1 || f_new != f || L₁_new != L₁ || L₂_new != L₂ || L₃_new != L₃ || κ₅_new != κ₅ || g_new != g || ν_new != ν
        println("ERROR: encountered simulation of multiple simultaneous temperatures or changed parameters.")
        exit()
    end

    if M_new != M
        println("Warning: one of the measurements has M = $(M_new) vs the previous M = $(M)")
    end

    jv = jackVars(energies -> specificHeat(energies, 1/temps[1]), E_by_T[1], num_blocks; skip_check=true)
    Cv_avg_by_T[i_f], var = jackEstimate(jv)
    Cv_err_by_T[i_f] = √(var)
    T_list[i_f] = temps[1]

    # Calculating average energy

    jv = jackVars(mean, E_by_T[1]./N, num_blocks; skip_check=true)
    E_avg_by_T[i_f], var = jackEstimate(jv)
    E_err_by_T[i_f] = √(var)

    # Calculating S⁺
    
    vo_di = JLD.load("vorticity.jld")
    S⁺_by_T = vo_di["sp"];
    # Measured k-vector components3
    kx = [two_pi/L₁*(x-1-L₁/2) for x = 1:L₁]
    ky = [two_pi/L₂*(y-1-L₂/2) for y = 1:L₂];
    # Normalizing structure function
    normalization = (L₁*L₂*f*two_pi)^2#(L^2*f*two_pi)^2
    S⁺_by_T = [S⁺_measures./normalization for S⁺_measures in S⁺_by_T];

    # Singeling out measurement series for a specific point (n₁, n₂)
    (n₁, n₂) = (Int(L₁/2)-1, 5+Int(L₂/2))
    hex_series = [fourier_lattice[n₁, n₂] for fourier_lattice in S⁺_by_T[1]]
    (n₁, n₂) = (Int(L₁/2)-3, 4+Int(L₂/2))
    quad_series = [fourier_lattice[n₁, n₂] for fourier_lattice in S⁺_by_T[1]]
    (n₁, n₂) = (Int(L₁/2)+3, 4+Int(L₂/2))
    stripe_series = [fourier_lattice[n₁, n₂] for fourier_lattice in S⁺_by_T[1]]
    (n₁, n₂) = (Int(L₁/2), 1+Int(L₂/2))
    control_series = [fourier_lattice[n₁, n₂] for fourier_lattice in S⁺_by_T[1]]

    jv_hex = jackVars(mean, hex_series, num_blocks; skip_check=true)
    jv_quad = jackVars(mean, quad_series, num_blocks; skip_check=true)
    jv_stripe = jackVars(mean, stripe_series, num_blocks; skip_check=true)
    jv_control = jackVars(mean, control_series, num_blocks; skip_check=true)

    hex_point_T[i_f], hex_var = jackEstimate(jv_hex)
    quad_point_T[i_f], quad_var = jackEstimate(jv_quad)
    stripe_point_T[i_f], stripe_var = jackEstimate(jv_stripe)
    control_T[i_f], control_var = jackEstimate(jv_control)

    hex_err_T[i_f] = √(hex_var)
    quad_err_T[i_f] = √(quad_var)
    stripe_err_T[i_f] = √(stripe_var)
    control_err_T[i_f] = √(control_var)

    # Calculating average amplitudes

    amplitude_di = JLD.load("amplitudes.jld");
    u⁻_avg_lattices = amplitude_di["um_xy"];
    u⁺_avg_lattices = amplitude_di["up_xy"];
    u⁺_by_T[i_f] = mean(u⁺_avg_lattices[1]./M_new)
    u⁻_by_T[i_f] = mean(u⁻_avg_lattices[1]./M_new)


    cd("../")
end

mkcd(plot_folder)


# Plotting amplitudes

plt = scatter(T_list, u⁺_by_T; ylabel=L"\langle\bar{u}^\pm\rangle", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label=L"u^+", xlims=(0,maximum(T_list)*1.05))
scatter!(plt, T_list, u⁻_by_T; label=L"u^-")
savefig(plt, "amplitudes_by_T.pdf")
plt = scatter(1 ./T_list, u⁺_by_T; ylabel=L"\langle\bar{u}^\pm\rangle", xlabel=L"\beta",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label=L"u^+", xlims=(0, maximum(1 ./T_list)*1.05))
scatter!(plt, 1 ./T_list, u⁻_by_T; label=L"u^-")
savefig(plt, "amplitudes_by_beta.pdf")

plt = scatter(T_list, abs.(u⁺_by_T.-u⁻_by_T); ylabel=L"\|\langle\bar{u}^+ - \bar{u}^-\rangle\|", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label="", xlims=(0,maximum(T_list)*1.05))
savefig(plt, "ampdiff_by_T.pdf")
plt = scatter(1 ./T_list, abs.(u⁺_by_T.-u⁻_by_T); ylabel=L"\|\langle\bar{u}^+ - \bar{u}^-\rangle\|", xlabel=L"\beta",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label="", xlims=(0, maximum(1 ./T_list)*1.05))
savefig(plt, "ampdiff_by_beta.pdf")

# Plotting structure factor

plt = scatter(T_list, stripe_point_T; yerror=stripe_err_T, ylabel=L"S^+(K)", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label=L"K_\mathrm{stripe} = \frac{2\pi}{L}[3, 3]", xlims=(0,maximum(T_list)*1.05))
scatter!(plt, T_list, hex_point_T; yerror=hex_err_T, label=L"K_\mathrm{hex} = \frac{2\pi}{L}[-1, 4]")
scatter!(plt, T_list, quad_point_T; yerror=quad_err_T, label=L"K_\mathrm{quad} = \frac{2\pi}{L}[-3, 3]")
savefig(plt, "stripe_by_T.pdf")

# Separate plot for only hex and quad points

plt = scatter(T_list, hex_point_T; yerror=hex_err_T, ylabel=L"S^+(K)", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label=L"K_\mathrm{hex} = \frac{2\pi}{L}[-1, 4]", xlims=(0,maximum(T_list)*1.05))
scatter!(plt, T_list, quad_point_T; yerror=quad_err_T, label=L"K_\mathrm{quad} = \frac{2\pi}{L}[-3, 3]")
savefig(plt, "hex_quad_by_T.pdf")
plt = scatter(1 ./T_list, hex_point_T; yerror=hex_err_T, ylabel=L"S^+(K)", xlabel=L"\beta",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              label=L"K_\mathrm{hex} = \frac{2\pi}{L}[-1, 4]", xlims=(0, maximum(1 ./T_list)*1.05))
scatter!(plt, 1 ./T_list, quad_point_T; yerror=quad_err_T, label=L"K_\mathrm{quad} = \frac{2\pi}{L}[-3, 3]")
savefig(plt, "hex_quad_by_beta.pdf")


# Plotting Cv

plt = scatter(T_list, Cv_avg_by_T./N, yerror=Cv_err_by_T./N; ylabel=L"\frac{C_v}{N}", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              ylims=(0, maximum(Cv_avg_by_T./N)*1.1), xlims=(0,maximum(T_list)*1.05), label="")
savefig(plt, "cv_by_T.pdf")
plt = scatter(1 ./ T_list, Cv_avg_by_T./N, yerror=Cv_err_by_T./N; ylabel=L"\frac{C_v}{N}", xlabel=L"\beta",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              ylims=(0, maximum(Cv_avg_by_T./N)*1.1), xlims=(0, maximum(1 ./T_list)*1.05), label="")
savefig(plt, "cv_by_beta.pdf")


# Plotting energy

plt = scatter(T_list, E_avg_by_T, yerror=E_err_by_T; ylabel=L"\langle E\rangle", xlabel="T",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              ylims=(0, maximum(E_avg_by_T)*1.1), xlims=(0,maximum(T_list)*1.05), label="")
savefig(plt, "E_by_T.pdf")
plt = scatter(1 ./T_list, E_avg_by_T, yerror=E_err_by_T; ylabel=L"\langle E\rangle", xlabel=L"\beta",
              title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)",
              ylims=(0, maximum(E_avg_by_T)*1.1), xlims=(0, maximum(1 ./T_list)*1.05), label="")
savefig(plt, "E_by_beta.pdf")


cd("../")
