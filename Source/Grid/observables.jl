############################################################################################################################
#                               Functions for observables
#__________________________________________________________________________________________________________________________#
############################################################################################################################
# Needs translate2DMat from utilities.jl for structureFunction
using LinearAlgebra # For dot(x,y)
using FFTW # For fast Fourier transform in structureFunction

# -------------------------------------------------------------------------------------------------
# Calculate gauge stiffness more effectively as a triplet
function dualStiffness(k::Tuple{T,T,T}, plaq_matrix::Array{Tuple{T,T,T}, 3}) where T<:Real
    gs_xx = Complex(0.0); gs_yy = Complex(0.0); gs_zz = Complex(0.0)
    L₁ = size(plaq_matrix, 1); L₂ = size(plaq_matrix, 2); L₃ = size(plaq_matrix, 3)

    for x = 1:L₁, y = 1:L₂, z = 1:L₃
        r = [x-1, y-1, z-1] #Spør 9cco om dette, tror det blir riktig
        gs_xx += plaq_matrix[x,y,z][1]*exp(im*dot(k,r))
        gs_yy += plaq_matrix[x,y,z][2]*exp(im*dot(k,r))
        gs_zz += plaq_matrix[x,y,z][3]*exp(im*dot(k,r))
    end

    norm = two_pi^2*L₁*L₂*L₃
    return abs2(gs_xx)/norm, abs2(gs_yy)/norm, abs2(gs_zz)/norm
end
function dualStiffnesses(cub::Cuboid)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    k₁ = (two_pi/L₁, 0.0, 0.0)
    k₂ = (0.0, two_pi/L₂, 0.0)
    k₃ = (0.0, 0.0, two_pi/L₃)
    
    plaq_matrix = fluxDensity(cub)
    ρˣˣₖ₂, _, ρᶻᶻₖ₂ = dualStiffness(k₂, plaq_matrix)
    ρˣˣₖ₃, ρʸʸₖ₃, _ = dualStiffness(k₃, plaq_matrix)
    _, ρʸʸₖ₁, ρᶻᶻₖ₁ = dualStiffness(k₁, plaq_matrix)
    
    ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₁, ρʸʸₖ₃, ρᶻᶻₖ₁, ρᶻᶻₖ₂
end

# -------------------------------------------------------------------------------------------------
# Calculates the specific heat using the temperature and a series of energies measured at this temp.
function specificHeat(energies::Array{T, 1}, β::T) where T<:Real
    return var(energies)*β^2
end

# -------------------------------------------------------------------------------------------------
# Calculates the Fourier Transform of 2D vortex lattices V⁺ and V⁻
# Uses fast fourier transform algorithm in the FFTW package and then translates the matrices
# so that we get 1st Brillouin zone. Assumes square matrices of equal size. Returns the same
# as calculating structureFunction (above) for all k = [k_x, k_y] where
# k_x, k_y ∈ 2π/L×[-L/2, L/2-1] and putting all the results in a matrix.
function structureFunction(V⁺::Array{R,2}, V⁻::Array{R,2}) where R<:Real
    L = size(V⁺, 1)
    S⁺_new = abs2.(bfft(V⁺)); S⁻_new = abs2.(bfft(V⁻))
    S⁺_new = translate2DMat(S⁺_new, Int(L/2)-1, Int(L/2)); S⁻_new = translate2DMat(S⁻_new, Int(L/2)-1, Int(L/2))
    S⁺_new, S⁻_new
end

