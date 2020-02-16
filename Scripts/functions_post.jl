# This script takes as input JLD files in the calling directory that starts with "energy_table".
# It is assumed these files contain a dictionary with an "Ens_by_T" and "T_list" entry that
# contains an Array{Array{Float64, 1}} of energies and Array{Float64, 1} of corresponding temperatures.
#
# Using this data, the script calculates the average specific heat with jackknife error,
# assuming a system size of N and using num_blocks blocks.
#
# It then uses the data for temperatures between T_min and T_max to reweigh the specific heat to the
# temperatures set in rw_Ts. All raw and reweighted cvs are saved to JLD data file.

using Distributions
using StatsBase
using LinearAlgebra
using LaTeXStrings
using BenchmarkTools
using Test
using Dates
using Primes
using MCMCDiagnostics

src_path = "/home/nicolai/mc/Source/Grid/"
include(src_path*"utilities.jl")
include(src_path*"plot_functions.jl")
ferr_path = src_path*"../FS/"
push!(LOAD_PATH, ferr_path)
using FerrenbergSwendsenReweighting

src_path = "/home/nicolai/mc/Source/"
include(src_path*"jackknife_estimates.jl")

#using Plots
#pyplot()

using JLD



##############################################################################################################
#
#                                   Cv Reweight Functions
#
##############################################################################################################
# Functions develped while trying to reweight Cv
# ____________________________________________________________________________________________________________

function specificHeat(energies::Array{T, 1}, β::T) where T<:Real
    return var(energies)*β^2
end

function loadEnergyFiles(; regex=r"energy_table.*\.jld")
    energy_files = filter(s -> isfile(s) && occursin(regex, s), readdir())
    println("Loading energy from $(length(energy_files)) files:")
    for file in energy_files; println(file); end
    N_F = length(energy_files)

    # Load temperatures and energy into memory
    # Energy and temperature storage
    E_pr_file = Array{Array{Array{Float64, 1}, 1}, 1}(undef, N_F)
    T_pr_file = Array{Array{Float64, 1}, 1}(undef, N_F)

    N_T = 0
    for (i_f, file) = enumerate(energy_files)
        temp_dict = JLD.load(file)
        E_pr_file[i_f] = temp_dict["Ens_by_T"]
        T_list = temp_dict["T_list"]
        N_T += length(T_list)
        T_pr_file[i_f] = T_list
    end

    # Then we flatten the vectors
    E_pr_temp = Array{Array{Float64, 1}, 1}(undef, N_T)
    T_list = Array{Float64, 1}(undef, N_T)
    j = 1
    for f = 1:N_F
        for i = 1:length(T_pr_file[f])
            E_pr_temp[j] = E_pr_file[f][i]
            T_list[j] = T_pr_file[f][i]
            j += 1
        end
    end
    E_pr_temp, T_list
end

function calculateCvUncertainty(es::Array{Array{R, 1}}, βs::Array{R, 1}, N::I; num_blocks=2^7) where {R<:Real, I<:Int}
    N_T = length(es)
    cv_avgs = Array{Float64, 1}(undef, N_T)
    cv_errs = Array{Float64, 1}(undef, N_T)
    for k = 1:N_T
        jv = jackVars(energies -> specificHeat(energies, βs[k])/N, es[k], num_blocks; skip_check=true)
        cv_avgs[k], var = jackEstimate(jv)
        cv_errs[k] = √(var)
    end
    cv_avgs, cv_errs
end

function restrictToInterval(Es::Array{Array{R, 1}, 1}, Ts::Array{R, 1}, T_lims::Tuple{R,R}) where R <: Real
    T_min, T_max = T_lims
    # First sorting lists
    perm = sortperm(Ts)
    Ts = Ts[perm]
    Es = Es[perm]

    res_indices = findall(T -> T>=T_min && T<=T_max, Ts)
    res_Ts = Ts[res_indices]
    res_Es = Es[res_indices]

    res_Es, res_Ts
end

function reweightCv(rw_βs::Array{R, 1}, Es::Array{Array{R, 1}, 1}, Ts::Array{R, 1}, N::I; rw = ReweightObj([1/T for T in Ts], Es; logarithms=true, verbose=true)) where {R <: Real, I <: Int}
    E²s = [[e^2 for e in e_list] for e_list in Es]
    rw_Es = reweight(Es, rw, rw_βs)
    rw_E²s = reweight(E²s, rw, rw_βs)
    rw_cvs = Array{R, 1}(undef, length(rw_Es))
    for i = 1:length(rw_Es)
        rw_cvs[i] = rw_βs[i]^2*(rw_E²s[i] - rw_Es[i]^2)/N
    end

    rw_cvs
end

function reweightCvWithUncertainty(rw_Ts, Ts, Es, N; num_blocks=2^7)
    βs = [1/T for T in Ts]
    rw_βs = [1/T for T in rw_Ts]

    # Setup reweighting object
    rw = ReweightObj(βs, Es; logarithms=true, verbose=true)
    init_guess = rw.ΔlogZs

    # This function is going to take in raw energies by T and give out reweighted cvs according to rw_βs
    function θ_estimator(os::Array{Array{T, 1}, 1}, es::Array{Array{R, 1}, 1}) where {T, R<:Real}
        # Create reweight object and solve FS equations
        rw = ReweightObj(βs, es; logarithms=true, initial_guess=init_guess, verbose=false);
        # Update initial guess.
        init_guess = rw.ΔlogZs

        reweightCv(rw_βs, es, Ts, N; rw=rw)
    end

    N_rw_T = length(rw_Ts)
    j_vars = jackVars(θ_estimator, Es, Es; num_blocks=num_blocks, skip_check=true)
    rw_Cvs = Array{Float64, 1}(undef, N_rw_T)
    rw_Cv_errs = Array{Float64, 1}(undef, N_rw_T)
    # Calculate jackknife variables using the reweighting estimator
    for k = 1:N_rw_T
        rw_Cvs[k], j_var = jackEstimate([j_vars[m][k] for m = 1:num_blocks])
        rw_Cv_errs[k] = √(j_var)
    end

    rw_Cvs, rw_Cv_errs
end

function findReweightedCvMax(rw_Ts::Array{R, 1}, Es::Array{Array{R, 1}}, Ts::Array{R, 1}, N::I; num_blocks=2^7, rw = ReweightObj([1/T for T in Ts], Es; logarithms=true, verbose=true)) where {R <: Real, I <: Int}
    βs = [1/T for T in Ts]
    rw_βs = [1/T for T in rw_Ts]
    init_guess = rw.ΔlogZs

    # This function is going to take in raw energies by T and give out T at index of max reweighted cv
    function θ_estimator(os::Array{Array{T, 1}, 1}, es::Array{Array{R, 1}, 1}) where {T, R<:Real}
        # Create reweight object and solve FS equations
        rw = ReweightObj(βs, es; logarithms=true, initial_guess=init_guess, verbose=false);
        # Update initial guess.
        init_guess = rw.ΔlogZs
        # Find reweighted cvs
        rw_cvs = reweightCv(rw_βs, es, Ts, N; rw=rw)
        # Return temperature at maximum reweighted Cv
        max_ind = findmax(rw_cvs)[2]
        rw_Ts[max_ind]
    end

    # Calculate jackknife variables using the reweighting estimator
    N_rw_T = length(rw_Ts)
    j_vars = jackVars(θ_estimator, Es, Es; num_blocks=num_blocks, skip_check=true)
    Tc, Tc_var = jackEstimate(j_vars)
    Tc_err = √(Tc_var)

    Tc, Tc_var, rw
end


function deeplength(A::AbstractArray)
    if typeof(A[1]) <: AbstractArray
        return sum([deeplength(x) for x in A])
    else
        return length(A)
    end
end

function reweightAndErrorCvForInterval(rw_Ts::Array{R, 1}, raw_Ts::Array{R, 1}, raw_Es::Array{Array{R, 1}, 1}, T_lims::Tuple{R, R}, N::I; num_blocks=2^7) where {R<:Real, I<:Int}
    T_min, T_max = T_lims
    # For reweighting we restrict the input-temperatures and energies to be between T_min and T_max.
    Es, Ts = restrictToInterval(raw_Es, raw_Ts, T_lims)
    N_T = length(Ts)
    βs = [1/T for T in Ts]

    # Calculate un-reweighted cv plot
    println("Calculating raw plot-values in reweighting interval.")
    cv_avgs, cv_errs = calculateCvUncertainty(Es, βs, N; num_blocks=num_blocks)

    println("Reweighting using $(N_T) temperature series with $(deeplength(Es)) energy measurements for $(length(rw_Ts)) reweighting temperatures.")
    rw_cvs, rw_cv_errs = reweightCvWithUncertainty(rw_Ts, Ts, Es, N; num_blocks=num_blocks)

    Ts, cv_avgs, cv_errs, rw_cvs, rw_cv_errs
end



##############################################################################################################
#
#                                   Post Temp Series Functions
#
##############################################################################################################
# Functions develped for plotting temperature series plots for cv, energy, S⁺, etc.
# ____________________________________________________________________________________________________________

function extractMetaSyst(;filename="meta.jld")

    meta_di = JLD.load(filename)
    m = meta_di["m"]; n = meta_di["n"]; L₃ = meta_di["L3"]; M = meta_di["M"]; κ₅ = meta_di["kap5"];
    L₁ = meta_di["L1"]; L₂ = meta_di["L2"]; g = meta_di["g"]; ν = meta_di["nu"]
    f = n/L₁ - m/L₂
    N = L₁*L₂*L₃;

    L₁, L₂, L₃, N, n, m, f, g, ν, κ₅
end
function extractMetaTemp(;filename="meta.jld")
    meta_di = JLD.load(filename)
    meta_di["T"]
end
function extractAmpMeasurements()
    amplitude_di = JLD.load("amplitudes.jld");
    # If we have old storage of up_avg where we forgot to specify up_avgs[k] then we have to figure out which
    # index we should use for the three possible
    u⁺_loads = amplitude_di["up_avg"]
    if typeof(u⁺_loads[1]) <: AbstractArray
        # Needed for old storage of u_avgs
        nus = [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        nu_dic = Dict([("$(ν)", mod(i-1,3)+1) for (i, ν) = enumerate(nus)])
        (_, _, _, _, _, _, _, _, ν, _)  = extractMetaSyst()
        i_ν = nu_dic["$(ν)"]
        u⁺s = u⁺_loads[i_ν]
        u⁻s = amplitude_di["um_avg"][i_ν]
        δu²s = amplitude_di["du2"][i_ν]
    else
        u⁺s = amplitude_di["up_avg"]
        u⁻s = amplitude_di["um_avg"]
        δu²s = amplitude_di["du2"]
    end
    u⁺s, u⁻s, δu²s
end
function searchExtendedPointsForMax(struct_lattices::Array{Array{R, 2}, 1}, point::Tuple{I,I}) where {R<:Real, I<:Int}
    x, y = point
    # First we expand the points searched through an extended number of points.
    extended_points = [(x, y), (mod(x+1-1,L₁)+1, y), (mod(x-1-1,L₁)+1, y), (x, mod(y+1-1,L₂)+1),
                       (mod(x+1-1,L₁)+1, mod(y+1-1,L₂)+1), (mod(x-1-1,L₁)+1, mod(y+1-1,L₂)+1),
                       (x, mod(y-1-1,L₂)+1), (mod(x+1-1,L₁)+1, mod(y-1-1,L₂)+1), (mod(x-1-1,L₁)+1, mod(y-1-1,L₂)+1)]
    extended_S⁺s = Array{Float64, 1}(undef, 9)
    # Calculate average at extended points
    for (e, e_point) = enumerate(extended_points)
        extended_S⁺s[e] = mean([fourier_lattice[e_point...] for fourier_lattice in struct_lattices])
    end
    # Only save the S⁺ series at the extended point with max average.
    # Return the index of the extended_point that gave the largest average.
    extended_points[findmax(extended_S⁺s)[2]]
end
function accumulateExtendedPoints(struct_lattices::Array{Array{R, 2}, 1}, point::Tuple{I,I}) where {R<:Real, I<:Int}
    x, y = point
    # First we expand the points searched through an extended number of points.
    extended_points = [(x, y), (mod(x+1-1,L₁)+1, y), (mod(x-1-1,L₁)+1, y), (x, mod(y+1-1,L₂)+1),
                       (mod(x+1-1,L₁)+1, mod(y+1-1,L₂)+1), (mod(x-1-1,L₁)+1, mod(y+1-1,L₂)+1),
                       (x, mod(y-1-1,L₂)+1), (mod(x+1-1,L₁)+1, mod(y-1-1,L₂)+1), (mod(x-1-1,L₁)+1, mod(y-1-1,L₂)+1)]
    [sum([lattice[e_point...] for e_point in extended_points]) for lattice in struct_lattices]
end
# Opens vorticity file, read the real space vorticity measurements, converts it into a separate measurement series
# for each point in "points"
function extractS⁺s(points::Array{Tuple{I,I}, 1}; fuzzy=false, accu=false) where I<:Int
    N_p = length(points)
    S⁺ss = Array{Array{Float64, 1}, 1}(undef, N_p)
    vo_di = JLD.load("vorticity.jld")
    V⁺_projs = vo_di["vp"]; V⁻_projs = vo_di["vm"]
    S⁺_meas = [structureFunction(V⁺_projs[m],V⁻_projs[m])[1] for m = 1:length(V⁺_projs)]
    L₁ = size(S⁺_meas[1], 1); L₂ = size(S⁺_meas[1], 2)

    for (p, point) = enumerate(S⁺_points)
        if accu
            S⁺ss[p] = accumulateExtendedPoints(S⁺_meas, point)
        else
            if fuzzy
                point = searchExtendedPointsForMax(S⁺_meas, point)
            end
            S⁺ss[p] = [fourier_lattice[point...] for fourier_lattice in S⁺_meas]
        end
    end
    S⁺ss
end

function extractS⁺_avg(points::Array{Tuple{I,I}, 1}) where I<:Int
    N_p = length(points)

    vo_di = JLD.load("vorticity.jld")
    S⁺_avg = vo_di["sp_avg"]
    [S⁺_avg[point...] for point in points]
end

function loadMeasurements(folders::Array{S, 1}, S⁺_points::Array{Tuple{I,I}, 1}; fuzzy=false, accu=false) where {I<:Int, S<:AbstractString}

    home_folder = pwd()
    N_F = length(folders)
    N_p = length(S⁺_points)

    Ess = Array{Array{Float64, 1}, 1}(undef, N_F)
    Ts = Array{Float64, 1}(undef, N_F)
    S⁺s_Ts_points = [Array{Array{Float64, 1}, 1}(undef, N_F) for p=1:N_p]
#    S⁺_avg_Ts_points = [Array{Float64, 1}(undef, N_F) for p=1:N_p]
    u⁺ss = Array{Array{Float64, 1}, 1}(undef, 0)
    u⁻ss = Array{Array{Float64, 1}, 1}(undef, 0)
    δu²ss = Array{Array{Float64, 1}, 1}(undef, 0)
    
    for (k, folder) = enumerate(folders)
        cd(folder)
        
        # Calculating specific heat

        ens_di = JLD.load("energies.jld")
        Ess[k] = ens_di["Es"]

        Ts[k] = extractMetaTemp()

        # Extracting S+ points

        S⁺ss = extractS⁺s(S⁺_points; fuzzy=fuzzy, accu=accu)
#        S⁺_avg_points = extractS⁺_avg(S⁺_points)
        for p = 1:N_p
            S⁺s_Ts_points[p][k] = S⁺ss[p]
#            S⁺_avg_Ts_points[p][k] = S⁺_avg_points[p]
        end

        # Extracting amplitudes

        if isfile("amplitudes.jld")
            u⁺s, u⁻s, δu²s = extractAmpMeasurements()
            push!(u⁺ss, u⁺s); push!(u⁻ss, u⁻s); push!(δu²ss, δu²s)
        end

        cd(home_folder)
    end

    Ts, Ess, S⁺s_Ts_points, u⁺ss, u⁻ss, δu²ss
end

function jackknifeMeanErr(os::Array{R, 1}; num_blocks=2^7, skip_check=true) where R<:Real
    jv = jackVars(mean, os, num_blocks; skip_check=skip_check)
    avg, var = jackEstimate(jv)
    err = √(var)
    avg, err
end

function jackknifeMeanErr(os::Array{Array{R, 1}, 1}; num_blocks=2^7) where R <: Real
    N_T = length(os)
    o_avgs = Array{R, 1}(undef, N_T)
    o_errs = Array{R, 1}(undef, N_T)

    for k = 1:N_T
        o_avgs[k], o_errs[k] = jackknifeMeanErr(os[k]; num_blocks=num_blocks)
    end
    o_avgs, o_errs
end

using FFTW

two_pi=2π
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

function getRelativeTuple(tuple::Tuple{I,I}, L₁::I, L₂::I) where I<:Int
    rel_tuple = tuple .- (Int(L₁/2), 1+Int(L₂/2))
    (rel_tuple[2], rel_tuple[1])
end

function pltTempSeries(Ts::Array{R, 1}, avgs::Array{R, 1}, errs::Array{R, 1}; ylabel::AbstractString="", x_margin = 0.005, ylims = (0, maximum(avgs)*1.1), title="", filename="obs_by_T.pdf") where R <: Real
    T_min = minimum(Ts); T_max = maximum(Ts)
    plt = scatter(Ts, avgs; yerror=errs, ylabel=ylabel, xlabel="T",
                  title=title, ylims=ylims, xlims=((1-x_margin)*T_min,T_max*(1+x_margin)), label="", xticks=T_min:0.1:T_max)
    savefig(plt, filename)
    nothing
end
function pltTempSeries(Ts::Array{R, 1}, avgs::Array{R, 1}; ylabel::AbstractString="", x_margin = 0.005, ylims = (0, maximum(avgs)*1.1), title="", filename="obs_by_T.pdf") where R <: Real
    T_min = minimum(Ts); T_max = maximum(Ts)
    plt = scatter(Ts, avgs; ylabel=ylabel, xlabel="T",
                  title=title, ylims=ylims, xlims=((1-x_margin)*T_min,T_max*(1+x_margin)), label="", xticks=T_min:0.1:T_max)
    savefig(plt, filename)
    nothing
end
function pltTempSeries(Ts::Array{R, 1}, avgss::Array{Array{R, 1}, 1}, errss::Array{Array{R, 1}, 1}; ylabel::AbstractString="", labels::Array{S, 1} = Array{S, 1}(undef, 0), x_margin = 0.005, y_margin = 0.1, title="", filename="obs_by_T.pdf") where {R <: Real, S <: AbstractString}
    T_min = minimum(Ts); T_max = maximum(Ts)
    plt = scatter(Ts, avgss[1]; yerror=errss[1], ylabel=ylabel, xlabel="T",
                  title=title, xlims=((1-x_margin)*T_min,T_max*(1+x_margin)), label=labels[1])
    if length(avgss) >= 2
        for i = 2:length(avgss)
            scatter!(plt, Ts, avgss[i]; yerror=errss[i], label=labels[i])
        end
    end
    savefig(plt, filename)
    nothing
end

