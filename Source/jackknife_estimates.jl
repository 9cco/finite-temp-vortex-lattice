# Types and function for calculating jackknife estimates
# Requires: divisors(::Int64) and mkcd(::AbstractString) from functions_msc.jl and the FerrenbergSwendsenReweighting package.
#
using Statistics
using Plots
gr()

####################################################################################################
#                            Types needed
#
####################################################################################################

mutable struct BlockInd
    N::Int64    # The number of elements in the time-series
    Nb::Int64   # The number of elements in a block
    n::Int64    # The number of blocks n = N/Nb
end

####################################################################################################
#                            Functions autocorrelation in blocks
#
####################################################################################################

function BlockInd(N::Int64, Nb::Int64)
    Nb != 0 || throw(error("Can not have 0 block length"))
    if N < Nb || N%Nb != 0
        throw(error("Domain error. $(N) not divisible by $(Nb)"))
    end
    n = N/Nb
    return BlockInd(N, Nb, n)
end

function iterator(B::BlockInd, t::Int64)
    return (t-1)*B.Nb+1:t*B.Nb
end

function mapToBlockedSeries(O::Array{T, 1}, B::BlockInd) where T<:Real
    O_aux = Array{T}(undef, B.n)
    for t = 1:B.n
        O_aux[t] = mean(O[iterator(B,t)])#(t-1)*B.Nb+1:t*B.Nb])
    end
    return O_aux
end

# Variance error based on the variance of the sample variance found on wikipedia. Basically we have
# err = √Var(sₓ²) = √( Var(s²/n) ) = √( 1/n²⋅Var(s²) ) = 1/n⋅√Var(s²), where s² is the unbiased sample variance
# estimator and sₓ² is the estimator of the variance of the sample mean.
# Then we use the formula for Var(s²) = 1/n⋅[μ₄ - (n-3)/(n-1)⋅s⁴]. 
function varianceError(Os::Vector{<:Real})
    n = length(Os)
    mn = sum(Os)/n
    # calculating estimator for fourth moment
    μ₄ = sum([(Oᵢ-mn)^4 for Oᵢ in Os])/n
    # calculating estimator for sample variance
    s² = var(Os)

    √((μ₄ - (n-3)/(n-1)*s²^2)/n)/n
end

function varByBlockLength(O::Array{T, 1}, Nb_list::Array{I, 1}) where {T<:Real, I<:Int}
    N = length(O)
    M = length(Nb_list)
    var_list = Array{Float64}(undef, M)
    var_error_list = Array{Float64}(undef, M)
    
    for (i, Nb) in enumerate(Nb_list)
        B = BlockInd(N,Nb)
        O_aux = mapToBlockedSeries(O, B)
        var_list[i]  = var(O_aux)/B.n
        var_error_list[i] = max(varianceError(O_aux), var_list[i]*√(2/(B.n-1)))
    end
    return var_list, var_error_list
end

# Checks if the element Y[i] is within the error of the points Y[j] for i<j such that 
# Y[i] ∈ [Y[j]-err[j], Y[j]+err[j]]
function isWithinErrorOfRightPoints(i::Int64, Y::Array{T,1}, err::Array{T,1}; err_weight=1.0) where T<:Real
    for j = i+1:length(Y)
        if Y[i] < Y[j]-err_weight*err[j] || Y[i] > Y[j]+err_weight*err[j]
            return false
        end
    end
    return true
end

function firstConvergedIndex(Y::Array{T,1}, err::Array{T,1}) where T<:Real
    for i = 1:length(Y)-1
        if isWithinErrorOfRightPoints(i, Y, err)
            return i
        end
    end
    return -1
end

function smallestBlockSize(O_series::Array{T,1}; write_plot=false, filename="block_size_variance_plt.pdf") where T<:Real
    N = length(O_series)
    # Find the possible block sizes of the series {O}
    # Remove last possible block size to make sure number of blocks > 2
    Nb_list = reverse(divisors(N))[1:end-1]
    
    # Calculate the variance its error when using the different possible block lengths
    var_list, error_list = varByBlockLength(O_series, Nb_list)
    if write_plot
        plt = scatter(Nb_list, var_list; xscale=:log2, yerror=error_list, yaxis="S²", xaxis="Nb",
            title="Variance when dividing the series using block length Nb",
            label="Total length ∼ 2^($(round(log2(N); digits=1)))")
        savefig(plt, filename)
    end

    i = firstConvergedIndex(var_list, error_list)
    # If a converged index could not be found it might mean that the autocorrelation time is comparable or larger than
    # the O_series length, thus we return a length larger than this.
    if i == -1
        return N+1
    end
    
    return Nb_list[i]
end

# Plots the variance plot of a single series of observables that could be correlated. In theory, the estimate of
# the variance should converge from below when the block length is much larger than the autocorrelation time
# of the time-series. num_blocks is the number of blocks intended to be used when dividing the series for jackknife
# analysis. Saves the plots to the current directory.
function singleSeriesVariancePlot(O_series::Vector{<:Real}, num_blocks::Int; 
                                  filename="block_size_variance_plt.pdf", title="Variance of blocks with length Nb")

    # Find the possible block sizes of the series {O}
    # Remove last possible block size to make sure number of blocks > 2
    N = length(O_series)
    Nb_list = reverse(divisors(N))[1:end-1]
    Nb_proposal = floor(Int64, N/num_blocks)
    extra = 1/3

    # Calculate the variance and its error when using the different possible block lengths
    var_list, error_list = varByBlockLength(O_series, Nb_list)
    mx = maximum(var_list)
    mi = minimum(var_list)
    δy = mx - mi
    plt = scatter(Nb_list, var_list; xscale=:log2, yerror=error_list, yaxis="S²", xaxis="Nb",
                  title=title,
                  label="Total length ∼ 2^($(round(log2(N); digits=1)))",
                  ylims = (mi-extra*δy, mx+extra*δy))
    # Add proposal line
    plot!(plt, [Nb_proposal]; xscale=:log2, seriestype=:vline, label="Proposed length ~ 2^($(round(log2(Nb_proposal); digits=1)))")
    savefig(plt, filename)
    nothing
end
# Plots variance vs block-length plots to current folder for multiple time-series in O_set using the function above.
function multipleSeriesVariancePlots(O_set::Vector{Vector{R}}, num_blocks::Int; file_ending="png",
                                     titles=["Variance of blocks with length Nb" for k = 1:length(O_set)]) where R<:Real

    num_series = length(O_set)
    for (k, O_series) = enumerate(O_set)
        singleSeriesVariancePlot(O_series, num_blocks; filename="variance_of_series_$(k)."*file_ending, title=titles[k])
    end
    nothing
end


####################################################################################################
#                            Functions the jackknife estimates proper
#
####################################################################################################

function low(B::BlockInd, t::Int64)
    return 1:(t-1)*B.Nb
end
function high(B::BlockInd, t::Int64)
    return t*B.Nb+1:B.N
end

# Returns Jₜ = {1,...,N} \ Bₜ = {1, ..., (t-1)*Nb} ∪ {t*Nb+1, ..., N}
function jackSeries(O_list::Array{T,1}, B::BlockInd, t::Int64) where T<:Real
    O_aux = Array{T}(undef, B.N-B.Nb)
    i = 0
    for i = low(B,t)
        O_aux[i] = O_list[i]
    end
    for j = high(B,t)
        O_aux[j-B.Nb] = O_list[j]
    end
    O_aux
end

function jackSet(O_set::Array{Array{T,1},1}, J_set::Array{BlockInd,1}, t::Int64) where T<:Real
    N₀ = length(O_set)
    reduced_time_series_set = Array{Array{T,1}}(undef, N₀)
    for k = 1:N₀
        reduced_time_series_set[k] = jackSeries(O_set[k], J_set[k], t)
    end
    return reduced_time_series_set
end

# Checks that the block number gives block sizes that are sufficiently large for all time-series such that
# autocorrelation effects are avoided.
function blockNumIsSmallEnough(n::Int64, O_set::Array{Array{T,1},1}; verbose=false) where T<:Real
    n>0 || throw(error("ERROR: The block size must be larger than 1 (n set to $(n))"))
    num_series = length(O_set)
    for k = 1:num_series
        # First we check that n divides the time-series length
        Nₖ = length(O_set[k])
        Nₖ%n == 0 || throw(error("ERROR: Number of blocks $(n) does not 
divide number of elements $(Nₖ) in time series $(k)"))
        Nbₖ = Int(Nₖ/n)
        Nb_smallest = smallestBlockSize(O_set[k])
        if Nbₖ < Nb_smallest
            if verbose
                smallestBlockSize(O_set[k]; write_plot=true)
                println("Problem in time-series k = $(k) with length $(Nₖ). 
Smallest block size returned $(Nb_smallest), while n = $(n) number of blocks yields block size $(Nbₖ).\nWriting variance plot.")
            end
            return false
        end
    end
    return true
end

function blockNumIsSmallEnough(n::Int64, O_set::Array{T,1}; verbose=false) where T<:Real
    return blockNumIsSmallEnough(n, [O_set]; verbose=verbose)
end

# Return jackknife variables assuming a function θ(O₁_set, O₂_set) that takes two corresponding sets of time-series. Each corresponding
# series in each set is reduced by the t'th block to create jackknife variable number t. The output of θ can be an array, in which case
# jackVars returns an array of an array of jackknife variables.
function jackVars(θ::Function, O1_set::Array{Array{T,1}, 1}, O2_set::Array{Array{T2,1}}; num_blocks=2^7, skip_check=true) where {T<:Real, T2}

    if skip_check   # If we skipped the check we still need to check if num_blocks is a possible number of blocks.
        length.(O1_set) .% num_blocks == fill(0, length(O1_set)) || throw(error("ERROR: The number of blocks $(num_blocks) does not divide
all series of first argument O1_set"))
        length.(O2_set) .% num_blocks == fill(0, length(O2_set)) || throw(error("ERROR: The number of blocks $(num_blocks) does not divide
all series of second argument O2_set"))
    else
        if !blockNumIsSmallEnough(num_blocks, O1_set; verbose=true)
            throw(error("ERROR: The number of blocks: $(num_blocks) was too large to reduce 
autocorrelation effects of independent blocks for all series in first argument set O1_set"))
        elseif !blockNumIsSmallEnough(num_blocks, O2_set; verbose=true)
            throw(error("ERROR: The number of blocks: $(num_blocks) was too large to reduce 
autocorrelation effects of independent blocks for all series in second argument set O2_set"))
        end
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number num_blocks.
    # This also checks that all series in the sets have lengths that are divisible by num_blocks.
    
    # Checking that the dimensions of O1_set and O2_set are the same so that the same J_set can be used for both
    length.(O1_set) == length.(O2_set) || throw(error("ERROR: Observable 1 set does not have the same dimensions as
observable 2"))
    
    N₀ = length(O1_set)
    # For each time-series we create a jackknife object
    J_set = Array{BlockInd}(undef, N₀)
    for k = 1:N₀
        Nₖ = length(O1_set[k])
        J_set[k] = BlockInd(Nₖ, Int(Nₖ/num_blocks))
    end
    
    # Then we calculate the jackknife variables
    j_var = θ(jackSet(O1_set, J_set, 1), jackSet(O2_set, J_set, 1))
    J_vars = Array{typeof(j_var), 1}(undef, num_blocks)
    J_vars[1] = j_var
    for t = 2:num_blocks
        J_vars[t] = θ(jackSet(O1_set, J_set, t), jackSet(O2_set, J_set, t))
    end
    
    # Return the finished set of Jackknife variables
    return J_vars
end


# Efficiency note: It seems that most of the time goes to creating all the different jackSeries. This could be
# prevented by only using the original series and only calculate the indices needed for the different jackknife series.
# But this will restrict the estimator θ to take both the original O_list as well as a list of indices that
# it should use. For θ that depends on O_set, then it should instead also take a set of lists of indices.

function jackVars(θ::Function, O_set::Array{Array{T,1}, 1}; num_blocks=2^7, skip_check=false) where T<:Real
    N₀ = length(O_set)
    if !skip_check && !blockNumIsSmallEnough(num_blocks, O_set; verbose=true)
        throw(error("ERROR: The number of blocks: $(num_blocks) was too large to reduce 
autocorrelation effects of independent blocks for all sets in { O_list[k] }"))
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number n.
    # This also checks that all series in the sets have lengths that are divisible by n.
    
    # For each time-series we create a jackknife object
    J_set = Array{BlockInd}(undef, N₀)
    for k = 1:N₀
        Nₖ = length(O_set[k])
        J_set[k] = BlockInd(Nₖ, Int(Nₖ/num_blocks))
    end
    
    # Then we calculate the jackknife variables
    j_var = θ(jackSet(O_set, J_set, 1))
    J_vars = Array{typeof(j_var), 1}(undef, num_blocks)
    J_vars[1] = j_var
    for t = 2:num_blocks
        J_vars[t] = θ(jackSet(O_set, J_set, t))
    end
    
    # Return the finished set of Jackknife variables
    return J_vars
end

function jackVars(θ::Function, O_list::Array{T,1}, n::Int64; skip_check=false) where T<:Real
    if !skip_check && !blockNumIsSmallEnough(n, O_list; verbose=true)
        throw(error("ERROR: The number of blocks n = $(n) was too large to reduce 
autocorrelation effects of independent blocks for O_list"))
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number n.
    # This also checks that all series in the sets have lengths that are divisible by n.
    
    # For each time-series we create a jackknife object
    N = length(O_list)
    J = BlockInd(N, Int(N/n))
    
    # Then we calculate the jackknife variables
    J_vars = Array{Real}(undef, n)
    for t = 1:n
        J_vars[t] = θ(jackSeries(O_list, J, t))
    end
    
    # Return the finished set of Jackknife variables
    return J_vars
end

# For this function we assume that θ takes care of only using the proper values of O_set, which are
# determined by the input J_set and t index.
function fastJackVars(θ::Function, O_set::Array{Array{T,1}}, n::Int64; skip_check=false) where T<:Real
    N₀ = length(O_set)
    if !skip_check && !blockNumIsSmallEnough(n, O_set; verbose=true)
        thow(error("ERROR: The number of blocks n = $(n) was too large to reduce 
autocorrelation effects of independent blocks for all sets in { O_list[k] }"))
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number n.
    # This also checks that all series in the sets have lengths that are divisible by n.
    
    # For each time-series we create a jackknife object
    J_set = Array{BlockInd}(undef, N₀)
    for k = 1:N₀
        Nₖ = length(O_set[k])
        J_set[k] = BlockInd(Nₖ, Int(Nₖ/n))
    end
    
    # Then we calculate the jackknife variables
    J_vars = Array{Real}(undef, n)
    for t = 1:n
        J_vars[t] = θ(O_set, J_set, t)
    end
    
    # Return the finished set of Jackknife variables
#    return mean(J_vars), (n-1)^2/n*var(J_vars)
    return J_vars
end

# For an array of Jackknife variables, return mean and variance
function jackEstimate(J_vars::Array{T,1}) where T<:Real
    n = length(J_vars)
    return mean(J_vars), (n-1)^2/n*var(J_vars)
end


####################################################################################################
#                            Functions using jackknife in reweighting
#
####################################################################################################

# Return an array of indices of series in O_set where smallest block size is larger than the
# proposed Nb
function findProblemSeries(O_set::Vector{Vector{T}}; num_blocks=2^7) where {T<:Real}

    N₀ = length(O_set)
    Nₖ = length.(O_set)
    problem_indices = Array{Int64, 1}(undef, 0)
    for (k, Os) = enumerate(O_set)
        Nb_prop = Nₖ[k]/num_blocks
        Nb_min = smallestBlockSize(Os)
        if Nb_prop < Nb_min
            push!(problem_indices, k)
        end
    end
    problem_indices
end

# If plot_all is true, this function plots the variance-by-block-length plots for all series in O_by_T and
# E_by_T in two separate folders called obs_folder and en_folder. The kind of plots produced can be determined
# with the file_ending keyword. If plot_all is false, then the function only produces the plots for series
# that has convergent block_size larger than the one dictated by num_blocks (number of blocks).
#   The function deletes the obs_folder and en_folder if clear is true. It also includes the temperatures
# in the title of the plots.
function varianceByBlockLengthPlots(O_by_T::Vector{Vector{T}}, E_by_T::Vector{Vector{R}}, β_list::Vector{R}; 
                                  obs_folder="var_by_block_length_obs", en_folder="var_by_block_length_en",
                                  file_ending="png", num_blocks=2^7, plot_all=false, clear=false) where {T, R<:Real}

    if clear
        rm(obs_folder; recursive=true)
        rm(en_folder; recursive=true)
    end

    if plot_all
        titles = ["Variance of blocks with length Nb, T=$(round(1/β; digits=3))" for β ∈ β_list]
        # We start by creating the observable plots
        mkcd(obs_folder)
        multipleSeriesVariancePlots(O_by_T, num_blocks; file_ending=file_ending, titles=titles)
        cd("../")

        # Then energy plots
        mkcd(en_folder)
        multipleSeriesVariancePlots(prob_series, num_blocks; file_ending=file_ending, titles=titles)
        cd("../")

    else
        # First we find series of observables where convergence is not found.
        prob_indices = findProblemSeries(O_by_T; num_blocks=num_blocks)
        if length(prob_indices) >= 1
            prob_series = [O_by_T[k] for k ∈ prob_indices]
            titles = ["Variance of blocks with length Nb, T=$(round(1/β_list[k]; digits=3))" for k ∈ prob_indices]
            mkcd(obs_folder)
            multipleSeriesVariancePlots(prob_series, num_blocks; file_ending=file_ending, titles=titles)
            cd("../")
        end

        # Then do the same for energies
        prob_indices = findProblemSeries(E_by_T; num_blocks=num_blocks)
        if length(prob_indices) >= 1
            prob_series = [E_by_T[k] for k ∈ prob_indices]
            titles = ["Variance of blocks with length Nb, T=$(round(1/β_list[k]; digits=1))" for k ∈ prob_indices]
            mkcd(en_folder)
            multipleSeriesVariancePlots(prob_series, num_blocks; file_ending=file_ending, titles=titles)
            cd("../")
        end
    end
    nothing
end


# Using jackknife variables to estimate the error when reweighting using multi-histogram reweighting in the
# FerrenbergSwendsenReweighting package.
#     Takes in a vector of (time)-series of observables O_by_T where each series corresponds to an element of inverse
# temperatures β in orig_βs. The E_by_T argument is also a vector of time series, but is of the energy the state had
# when observable O_by_T[k][i] was measured, i.e. E_by_T[k][i] corresponds to O_by_T[k][i] and orig_βs[k] is the inverse
# temperature at which the time-series O_by_T[k] and E_by_T[k] were measured. The energies as then used to reweight
# the observables to a new (inverse) temperature range at rwt_βs. The error in the new reweighted observables is  found
# by creating jackknife blocks of the O_by_T and E_by_T sets.
function reweightWithError(O_by_T::Vector{Vector{T}}, E_by_T::Vector{Vector{R}}, orig_βs::Vector{R},
        rwt_βs::Vector{R}; num_blocks=2^7, skip_check=true, make_plots=true, file_ending="png") where {T, R<:Real}
    
    (N₀ = length(orig_βs)) == length(O_by_T) == length(E_by_T) || throw(error("ERROR: Length of input vectors did not agree.
Are you sure that we have corresponding observables, energies and inverse temperatures?"))
    length.(O_by_T) == length.(E_by_T) || throw(error("ERROR: Length of time series of observables 
and energies did not agree"))

    if make_plots    # Make plots of variance of non-convergent series, both for observable and energies.
        varianceByBlockLengthPlots(O_by_T, E_by_T, orig_βs; num_blocks=num_blocks, file_ending=file_ending)
    end

    N_temps = length(rwt_βs)
    rwt_Os = Array{T, 1}(undef, N_temps)
    rwt_errs = Array{T, 1}(undef, N_temps)
    
    # Create reweight object based on full set of measurements for the purpose of an initial guess of ΔlogZs
    rw = ReweightObj(orig_βs, E_by_T; logarithms=true, verbose=true);
    initial_guess = rw.ΔlogZs
    
    # Estimator function takes in observable array and corresponding energies,
    # giving reweighted estimates for each β in rwt_βs. Used in a jackknife estimate, this solves the
    # FS equations for a reduced set of energies for each jackknife variable.
    # Assumes that orig_βs, rwt_βs and initial_guess are defined in outer scope.
    function θ_estimator(os::Array{Array{T, 1}, 1}, es::Array{Array{R, 1}, 1}) where {T, R<:Real}
        N₀ = length(os)
        # Create reweight object and solve FS equations
        rw = ReweightObj(orig_βs, es; logarithms=true, initial_guess=initial_guess, verbose=false);
        # Update initial guess.
        initial_guess = rw.ΔlogZs
        # Return vector of reweighted observables.
        reweight(os, rw, rwt_βs)
    end
    
    # Calculate jackknife variables using the θ-estimator giving a num_blcoks length vector of
    # vectors, where each inner vector contains the jackknife variable for each β in rwt_βs.
    j_vars = jackVars(θ_estimator, O_by_T, E_by_T; num_blocks=num_blocks, skip_check=skip_check);
    
    # Now we use the jackknife-variables to estimate the mean and error in the measurements.
    ot = typeof(O_by_T[1][1])
    rwt_Os = Array{ot, 1}(undef, N_temps)
    rwt_errs = Array{ot, 1}(undef, N_temps)
    for i = 1:N_temps
        j_mean, j_var = jackEstimate([j_vars[m][i] for m = 1:num_blocks])
        rwt_Os[i] = j_mean
        rwt_errs[i] = √(j_var)
    end
    
    rwt_Os, rwt_errs
end
