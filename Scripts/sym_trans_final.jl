# We assume the script is called from a folder with a final state .jld file
# It takes the state and performes symmetry transformations like C₄ and C₄⁻¹,
# calculates the resulting energy and plots the increase in energy given
# different parts of the energy.

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

@everywhere src_path = "/home/nicolai/mc/Source/Grid/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModuleUpdate

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")
include(src_path*"plot_functions.jl")

@everywhere src_path = "/home/nicolai/mc/Source/"
@everywhere include(src_path*"jackknife_estimates.jl")

using Plots
pyplot()

using JLD

# Make and enter folder for storing results
plot_path = "SymmTrans"
mkcd(plot_path)
cd("../")

# Collecting and printing meta info

meta_di = JLD.load("meta.jld")
n = meta_di["n"]; m = meta_di["m"]; M = meta_di["M"]; T = meta_di["T"]; κ₅ = meta_di["kap5"]; κ = meta_di["kap"]
Δt = meta_di["dt"]; L₁ = meta_di["L1"]; L₂ = meta_di["L2"]; L₃ = meta_di["L3"]; g = meta_di["g"]; ν = meta_di["nu"]
M_amp = meta_di["M_amp"];
N_g = 1

f = n/L₁ - m/L₂
println("Calculating rotated energy differences of state with")
println("\nfL₁ = $(f*L₁), L₁ = $(L₁), L₂ = $(L₂), L₃ = $(L₃), g = $(g)")
println("T = $(T), ν = $(ν), κ₅ = $(κ₅), κ = $(κ), n = $(n), m = $(m)")
T_round = round(T; digits=2)
N = L₁*L₂*L₃;


init_file = "final_state_g=$(round(g; digits=3))_nu=$(ν)_kap=$(κ).jld"
control_syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,κ,n,m)
split = (1,1,1)
fixRC()
println("Constructing file from $(init_file)")

# Reading the first initial file.
init_file_di = JLD.load(init_file)
init_lattice = init_file_di["lattice"]
init_syst = init_file_di["syst"]
init_T = init_file_di["T"]
init_controls = init_file_di["controls"]

if init_syst != control_syst
    println("ERROR: Parameters in initial-states file does not match target parameters in script.")
    exit(-1)
end

function rotateXYCCW(A::Array{LatticeSite, 3})
    # First we rotated the lattice site positions
    rotated_A = copy.(permutedims(A, [2 1 3])[end:-1:1, :,:])
    # Then on each position we have to update the link variables.
    L₁ = size(rotated_A, 1); L₂ = size(rotated_A, 2); L₃ = size(rotated_A, 3)
    # To move the values of A around we create a temporary lattice for each A where
    # the values are stored correctly.
    rot_A₁s = Array{Float64, 3}(undef, L₁,L₂,L₃)
    rot_A₂s = Array{Float64, 3}(undef, L₁,L₂,L₃)
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rot_A₁s[x,y,z] = -rotated_A[mod(x-1+1,L₁)+1, y,z].A₂
        rot_A₂s[x,y,z] = rotated_A[x,y,z].A₁
    end
    # Then we store the corret A values into the lattice sites
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rotated_A[x,y,z].A₁ = rot_A₁s[x,y,z]
        rotated_A[x,y,z].A₂ = rot_A₂s[x,y,z]
    end
    # Then on each position we let θ₁ → θ₂ + π; θ₂ → θ₁
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        θ₁ = rotated_A[x,y,z].θ⁺
        rotated_A[x,y,z].θ⁺ = rotated_A[x,y,z].θ⁻ + π
        rotated_A[x,y,z].θ⁻ = θ₁
    end
    # Finally we need to update the position values
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rotated_A[x,y,z].x = x
        rotated_A[x,y,z].y = y
    end
    rotated_A
end
function rotateXYCCW(s::SystConstants)
    SystConstants(s.L₁, s.L₂, s.L₃, s.g⁻²,
    s.ν, s.κ₅, s.κ, -s.m, -s.n)
end
function rotateXYCW(A::Array{LatticeSite, 3})
    rotated_A = copy.(permutedims(A, [2 1 3])[:, end:-1:1, :])
    # Then on each position we have to update the link variables.
    L₁ = size(rotated_A, 1); L₂ = size(rotated_A, 2); L₃ = size(rotated_A, 3)
    # To move the values of A around we create a temporary lattice for each A where
    # the values are stored correctly.
    rot_A₁s = Array{Float64, 3}(undef, L₁,L₂,L₃)
    rot_A₂s = Array{Float64, 3}(undef, L₁,L₂,L₃)
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rot_A₁s[x,y,z] = rotated_A[x,y,z].A₂
        rot_A₂s[x,y,z] = -rotated_A[x, mod(y-1+1,L₂)+1,z].A₁
    end
    # Then we store the corret A values into the lattice sites
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rotated_A[x,y,z].A₁ = rot_A₁s[x,y,z]
        rotated_A[x,y,z].A₂ = rot_A₂s[x,y,z]
    end
    # Then on each position we let θ₁ → θ₂; θ₂ → θ₁ + π
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        θ₁ = rotated_A[x,y,z].θ⁺
        rotated_A[x,y,z].θ⁺ = rotated_A[x,y,z].θ⁻
        rotated_A[x,y,z].θ⁻ = θ₁ + π
    end
    # Finally we need to update the position values
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        rotated_A[x,y,z].x = x
        rotated_A[x,y,z].y = L₂ - y + 1
    end
    rotated_A
end
function rotateXYCW(s::SystConstants)
    SystConstants(s.L₁, s.L₂, s.L₃, s.g⁻²,
    s.ν, s.κ₅, s.κ, -s.m, s.n)
end
function transposeXY(A::Array{T, 3}) where T
    copy.(permutedims(A, [2 1 3]))
end

import CuboidModuleUpdate.SystConstants
function setupNeighbors(lattice::Array{LatticeSite, 3}, pos::Tuple{Int64,Int64,Int64}, L₁::Int64,L₂::Int64,L₃::Int64)
    x,y,z = pos
    
    # Setup neighbors with periodic boundary conditions
    ϕᵣ₊₁ = lattice[mod(x-1+1,L₁)+1, y, z]
    ϕᵣ₋₁ = lattice[mod(x-1-1,L₁)+1, y, z]
    ϕᵣ₊₂ = lattice[x, mod(y-1+1,L₂)+1, z]
    ϕᵣ₋₂ = lattice[x, mod(y-1-1,L₂)+1, z]
    ϕᵣ₊₃ = lattice[x, y, mod(z-1+1,L₃)+1]
    ϕᵣ₋₃ = lattice[x, y, mod(z-1-1,L₃)+1]
    NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
end 

# We assume density is defined on the form
# function density(::LatticeSite, ::NearestNeighbors, ::SystConstants)
function energy(lattice::Array{LatticeSite,3}, transformed_lattice::Array{LatticeSite,3},
        density::Function, syst::SystConstants, trans_syst::SystConstants)
    
    L₁ = syst.L₁; L₂ = syst.L₂; L₃ = syst.L₃
    @test L₁ == size(lattice, 1) && L₂ == size(lattice, 2) && L₃ == size(lattice, 3)
    energy = 0.0
    energy′ = 0.0
    
    # Loop through the lattice
    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        ϕ = lattice[x,y,z]
        ϕ′ = transformed_lattice[x,y,z]
        
        # Setup neighbors with periodic boundary conditions
        nb = setupNeighbors(lattice, (x,y,z), L₁,L₂,L₃)
        nb′ = setupNeighbors(transformed_lattice, (x,y,z), L₁,L₂,L₃)
        
        # Calculated energy given supplied density
        energy += density(ϕ, nb, syst)
        energy′ += density(ϕ′, nb′, trans_syst)
    end
    
    energy, energy′
end
import CuboidModuleUpdate.NearestNeighbors
import CuboidModuleUpdate.linkVariables

function fₖ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    ϕᵣ₊₁ = nb.ϕᵣ₊₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Complete kinetic term.
    Fₖ = 2*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁺-A₁)
          + (ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁺-A₂)
    + c.κ₅*((ϕ.u⁺)^2 - ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺-ϕ.θ⁺-A₃))
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁻-A₁)
          + (ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁻-A₂)
    + c.κ₅*((ϕ.u⁻)^2 - ϕ.u⁻*ϕᵣ₊₃.u⁻*cos(ϕᵣ₊₃.θ⁻-ϕ.θ⁻-A₃)))
    energy = Fₖ
end
function fᵥ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    ϕᵣ₊₁ = nb.ϕᵣ₊₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Potential energy term
    Fᵥ = (0.5*(1+(1+c.ν)/2)*(ϕ.u⁺^4 + ϕ.u⁻^4) 
          + 0.5*(1-c.ν)*(ϕ.u⁺*ϕ.u⁻)^2*(2+cos(2*(ϕ.θ⁺-ϕ.θ⁻))) 
          - (ϕ.u⁺^2 + ϕ.u⁻^2))
    energy = Fᵥ
end
function fₐₙ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    ϕᵣ₊₁ = nb.ϕᵣ₊₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Anisotropi terms
    Fₐₙ = (c.ν+1)*(ϕᵣ₊₁.u⁻*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - A₁)
                 - ϕᵣ₊₂.u⁻*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - A₂)
                 + ϕᵣ₊₂.u⁺*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - A₂)
                 - ϕᵣ₊₁.u⁺*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - A₁))

    energy = Fₐₙ
end
function fₘ(ϕ::LatticeSite,nb::NearestNeighbors, c::SystConstants)
    energy = 0.0
    A₁, A₂, A₃ = linkVariables(ϕ, c)
    ϕᵣ₊₁ = nb.ϕᵣ₊₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    # Mixed gradient terms
    Fₘ = c.κ*(1-c.ν)*(ϕᵣ₊₁.u⁺*ϕᵣ₊₂.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕᵣ₊₂.θ⁻-A₁+A₂)
                -ϕᵣ₊₁.u⁺*ϕ.u⁻*cos(ϕᵣ₊₁.θ⁺-ϕ.θ⁻-A₁)
                -ϕᵣ₊₂.u⁻*ϕ.u⁺*cos(ϕᵣ₊₂.θ⁻-ϕ.θ⁺-A₂)
                +ϕᵣ₊₂.u⁺*ϕᵣ₊₁.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕᵣ₊₁.θ⁻-A₂+A₁)
                -ϕᵣ₊₁.u⁻*ϕ.u⁺*cos(ϕᵣ₊₁.θ⁻-ϕ.θ⁺-A₁)
                -ϕᵣ₊₂.u⁺*ϕ.u⁻*cos(ϕᵣ₊₂.θ⁺-ϕ.θ⁻-A₂)
                +2*ϕ.u⁺*ϕ.u⁻*cos(ϕ.θ⁺-ϕ.θ⁻))
    energy = Fₘ
end
# Extending maxwell to work with the above machinery
function maxwell(ϕ::LatticeSite, nb::NearestNeighbors, syst::SystConstants)
    CuboidModuleUpdate.maxwell(ϕ, nb.ϕᵣ₊₁, nb.ϕᵣ₊₂, nb.ϕᵣ₊₃, syst)
end
function plotTransformationEnergies(lattice::Array{LatticeSite, 3}, trans_lattice::Array{LatticeSite, 3}, 
        syst::SystConstants, trans_syst::SystConstants;
        fr_title="Transformed relative energies", fr_name="transformed_fr.pdf", max_name="transformation_maxwell.pdf", 
        plot_path="SymmTrans")
    
    # Now we may calculate the energy differences between rotated and normal lattice for different energy densities.
    trans_percentages = Array{Float64, 1}(undef, 0)
    nul_tot_en, rot_en = energy(lattice, trans_lattice, CuboidModuleUpdate.fᵣ, syst, trans_syst)
    push!(trans_percentages, (rot_en-nul_tot_en)/nul_tot_en*100)
    nul_en, rot_en = energy(lattice, trans_lattice, fₖ, syst, trans_syst)
    push!(trans_percentages, (rot_en-nul_en)/nul_tot_en*100)
    nul_en, rot_en = energy(lattice, trans_lattice, fᵥ, syst, trans_syst)
    push!(trans_percentages, (rot_en-nul_en)/nul_tot_en*100)
    nul_en, rot_en = energy(lattice, trans_lattice, fₐₙ, syst, trans_syst)
    push!(trans_percentages, (rot_en-nul_en)/nul_tot_en*100)
    nul_en, rot_en = energy(lattice, trans_lattice, fₘ, syst, trans_syst)
    push!(trans_percentages, (rot_en-nul_en)/nul_tot_en*100)
    plt = bar(trans_percentages; xticks = ([1:5...], ["tot", "fₖ", "fᵥ", "fₐₙ", "fₘ"]),
        ylabel="(E_t - E)/E_total %", title=fr_title, label="")
    savefig(plt, plot_path*"/"*fr_name)
    
    # Calulating maxwell vs. en. density
    maxwell_percentages = Array{Float64, 1}(undef, 0)
    nul_en, rot_en = energy(lattice, trans_lattice, CuboidModuleUpdate.fᵣ, syst, trans_syst)
    nul_max, rot_max = energy(lattice, trans_lattice, maxwell, syst, trans_syst)
    push!(maxwell_percentages, (rot_en + rot_max - nul_en - nul_max)/(nul_en+nul_max)*100, 
        (rot_en-nul_en)/(nul_en+nul_max)*100, (rot_max-nul_max)/(nul_en+nul_max)*100)
    plt = bar(maxwell_percentages; xticks = ([1:3...], ["tot", "fᵣ", "maxwell"]),
        ylabel="(E_t - E)/E_total %", title=fr_title, label="")
    savefig(plt, plot_path*"/"*max_name)
    
    # Saving percentages
    open(plot_path*"/"*fr_title*".txt", "w") do io
        writedlm(io, [trans_percentages, maxwell_percentages])
    end
    
    trans_percentages, maxwell_percentages
end

# Now we use the above functions to rotate the lattice and calculate energy for plotting

# Now we may calculate the energy differences between rotated and normal lattice for different energy densities.
rotated_lattice = rotateXYCCW(init_lattice)
rot_syst = rotateXYCCW(init_syst)
plotTransformationEnergies(init_lattice, rotated_lattice, init_syst, rot_syst; fr_title="CCW. rotation",
    fr_name="ccw_rotation_fr.pdf", max_name="ccw_rotation_maxwell.pdf", plot_path=plot_path)

cw_rotated_lattice = rotateXYCW(init_lattice)
cw_rot_syst = rotateXYCW(init_syst)
plotTransformationEnergies(init_lattice, cw_rotated_lattice, init_syst, cw_rot_syst; fr_title="CW. rotation",
    fr_name="cw_rotation_fr.pdf", max_name="cw_rotation_maxwell.pdf", plot_path=plot_path)
