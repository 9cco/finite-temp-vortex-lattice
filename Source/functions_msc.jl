####################################################################################################
#                            Misc. help functions
#
####################################################################################################

# --------------------------------------------------------------------------------------------------
# Helping function for calculating position vector in [x,y] format from a lattice position
# assuming origo is in lower left corner of the square lattice. L is the size of the lattice along one dimension.
function getVectorPosition(L::Int64, pos::Array{Int64,1})
    v_pos = pos[1]
    h_pos = pos[2]
    return [h_pos-1, L-v_pos]
end

# --------------------------------------------------------------------------------------------------
function autocorrTime{T<:Real}(O::Array{T,1}, c::Float64=5.0)
    N = size(O,1)
    
    # Estimate all correlation functions for different lags=[0:N]
    ρ = autocor(O, 0:(N-1))
    # Estimate correlation times. Each element in the array is an estimate for the correlation time τ(m) 
    # where summing over the i=1:M lagged autocorrelation functions using different numbers  M <= N
    τ = 2.*cumsum(ρ).-1
    
    for m=1:N
        if m < c*τ[m]
            # We return the estimate for τ that is such that M < c*τ(m) 
            return τ[m]
        end
    end
    return τ[N]
end

# --------------------------------------------------------------------------------------------------
# Finds the highest divisor of M that is less than or equal to N.
# If none of the divisors of M are less than N, then 
# We assume that M,N>1
function highestDivisor(M::Int, N::Int)
    m = M
    if m <= N
        return m
    end
    # Now we know that M > N
    # First we calculate all the divisors of M except for 1 and M
    # and save them in the divisors array.
    pr = primes(2,Int(floor(sqrt(M)))+2)
    L = size(pr, 1)
    divisors = Array{Int}(0)
    for i = 1:L
        while(m%pr[i]==0)
            m /= pr[i]
            if m == 1
                break
            end
            push!(divisors, m)
            push!(divisors, M/m)
        end
    end
    
    # Then we sort this array so that high divisors come first and
    # iterate through it until a match is found.
    divisors = sort(divisors, rev=true)
    for div in divisors
        if div <= N
            return div
        end
    end
    
end
