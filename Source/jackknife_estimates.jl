# Types and function for calculating jackknife estimates
# Requires: divisors(::Int64) from functions_msc.jl
#

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

function mapToBlockedSeries{T<:Real}(O::Array{T, 1}, B::BlockInd)
    O_aux = Array{T}(B.n)
    for t = 1:B.n
        O_aux[t] = mean(O[iterator(B,t)])#(t-1)*B.Nb+1:t*B.Nb])
    end
    return O_aux
end

function varByBlockLength{T<:Real, I<:Int}(O::Array{T, 1}, Nb_list::Array{I, 1})
    N = length(O)
    M = length(Nb_list)
    var_list = Array{Float64}(M)
    var_error_list = Array{Float64}(M)
    
    for (i, Nb) in enumerate(Nb_list)
        B = BlockInd(N,Nb)
        O_aux = mapToBlockedSeries(O, B)
        var_list[i]  = var(O_aux)/B.n
        var_error_list[i] = var_list[i]*√(2/(B.n-1))
    end
    return var_list, var_error_list
end

# Checks if the element Y[i] is within the error of the points Y[j] for i<j such that 
# Y[i] ∈ [Y[j]-err[j], Y[j]+err[j]]
function isWithinErrorOfRightPoints{T<:Real}(i::Int64, Y::Array{T,1}, err::Array{T,1}; err_weight=1.0)
    for j = i+1:length(Y)
        if Y[i] < Y[j]-err_weight*err[j] || Y[i] > Y[j]+err_weight*err[j]
            return false
        end
    end
    return true
end

function firstConvergedIndex{T<:Real}(Y::Array{T,1}, err::Array{T,1})
    for i = 1:length(Y)-1
        if isWithinErrorOfRightPoints(i, Y, err)
            return i
        end
    end
    return -1
end

function smallestBlockSize{T<:Real}(O_series::Array{T,1}; write_plot=false, filename="block_size_variance_plt.png")
    N = length(O_series)
    # Find the possible block sizes of the series {O}
    Nb_list = reverse(divisors(N))
    
    # Calculate the variance its error when using the different possible block lengths
    var_list, error_list = varByBlockLength(O_series, Nb_list)
    i = firstConvergedIndex(var_list, error_list)
    if i == -1
        plt = scatter(Nb_list, var_list; xscale=:log2, yerror=error_list, yaxis="S²", xaxis="Nb",
            title="Variance when dividing the series using block length Nb")
        savefig(plt, filename)
        throw(error("ERROR: Could not determine minimum block size. Increase series length?"))
    end
    if write_plot
        plt = scatter(Nb_list, var_list; xscale=:log2, yerror=error_list, yaxis="S²", xaxis="Nb",
            title="Variance when dividing the series using block length Nb")
        savefig(plt, filename)
    end
    return Nb_list[i]
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
function jackSeries{T<:Real}(O_list::Array{T,1}, B::BlockInd, t::Int64)
    O_aux = Array{T}(B.N-B.Nb)
    i = 0
    for i = low(B,t)
        O_aux[i] = O_list[i]
    end
    for j = high(B,t)
        O_aux[j-B.Nb] = O_list[j]
    end
    O_aux
end

function jackSet{T<:Real}(O_set::Array{Array{T,1},1}, J_set::Array{BlockInd,1}, t::Int64)
    N₀ = length(O_set)
    reduced_time_series_set = Array{Array{T,1}}(N₀)
    for k = 1:N₀
        reduced_time_series_set[k] = jackSeries(O_set[k], J_set[k], t)
    end
    return reduced_time_series_set
end

# Checks that the block number gives block sizes that are sufficiently large for all time-series such that
# autocorrelation effects are avoided.
function blockNumIsSmallEnough{T<:Real}(n::Int64, O_set::Array{Array{T,1},1}; verbose=false)
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

function blockNumIsSmallEnough{T<:Real}(n::Int64, O_set::Array{T,1}; verbose=false)
    return blockNumIsSmallEnough(n, [O_set]; verbose=verbose)
end

# Efficiency note: It seems that most of the time goes to creating all the different jackSeries. This could be
# prevented by only using the original series and only calculate the indices needed for the different jackknife series.
# But this will restrict the estimator θ to take both the original O_list as well as a list of indices that
# it should use. For θ that depends on O_set, then it should instead also take a set of lists of indices.

function jackVars{T<:Real}(θ::Function, O_set::Array{Array{T,1}}, n::Int64; skip_check=false)
    N₀ = length(O_set)
    if !skip_check && !blockNumIsSmallEnough(n, O_set; verbose=true)
        throw(error("ERROR: The number of blocks n = $(n) was too large to reduce 
autocorrelation effects of independent blocks for all sets in { O_list[k] }"))
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number n.
    # This also checks that all series in the sets have lengths that are divisible by n.
    
    # For each time-series we create a jackknife object
    J_set = Array{BlockInd}(N₀)
    for k = 1:N₀
        Nₖ = length(O_set[k])
        J_set[k] = BlockInd(Nₖ, Int(Nₖ/n))
    end
    
    # Then we calculate the jackknife variables
    J_vars = Array{Real}(n)
    for t = 1:n
        J_vars[t] = θ(jackSet(O_set, J_set, t))
    end
    
    # Return the finished set of Jackknife variables
#    return mean(J_vars), (n-1)^2/n*var(J_vars)
    return J_vars
end

function jackVars{T<:Real}(θ::Function, O_list::Array{T,1}, n::Int64; skip_check=false)
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
    J_vars = Array{Real}(n)
    for t = 1:n
        J_vars[t] = θ(jackSeries(O_list, J, t))
    end
    
    # Return the finished set of Jackknife variables
    return J_vars
end

# For this function we assume that θ takes care of only using the proper values of O_set, which are
# determined by the input J_set and t index.
function fastJackVars{T<:Real}(θ::Function, O_set::Array{Array{T,1}}, n::Int64; skip_check=false)
    N₀ = length(O_set)
    if !skip_check && !blockNumIsSmallEnough(n, O_set; verbose=true)
        thow(error("ERROR: The number of blocks n = $(n) was too large to reduce 
autocorrelation effects of independent blocks for all sets in { O_list[k] }"))
    end
    # Now we have checked that the jackknife estimate will work for the chosen block number n.
    # This also checks that all series in the sets have lengths that are divisible by n.
    
    # For each time-series we create a jackknife object
    J_set = Array{BlockInd}(N₀)
    for k = 1:N₀
        Nₖ = length(O_set[k])
        J_set[k] = BlockInd(Nₖ, Int(Nₖ/n))
    end
    
    # Then we calculate the jackknife variables
    J_vars = Array{Real}(n)
    for t = 1:n
        J_vars[t] = θ(O_set, J_set, t)
    end
    
    # Return the finished set of Jackknife variables
#    return mean(J_vars), (n-1)^2/n*var(J_vars)
    return J_vars
end

# For an array of Jackknife variables, return mean and variance
function jackEstimate{T<:Real}(J_vars::Array{T,1})
    return mean(J_vars), (n-1)^2/n*var(J_vars)
end
