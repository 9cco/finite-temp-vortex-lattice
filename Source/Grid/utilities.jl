using StatsBase
using MCMCDiagnostics
using Primes

# --------------------------------------------------------------------------------------------------
# Generates a geometric series xᵢ where xᵢ = x₁ for i = 1 and xᵢ = xₛ for i = N. If the exponent
# is set to something other than 1, then xᵢ follows a generalized geometric series xᵢ ∼ (xₛ/x₁)^(i*exponent)
function geometricSeries(x₁::R, xₛ::R, N::Int; exponent=1) where R<:Real
    frac = (xₛ/x₁)
    f(i) = frac^(exponent*(i-1)/(N-1))*x₁*(1 - frac)/(1-frac^exponent) + xₛ*(1-frac^(exponent-1))/(1-frac^exponent)
    return [f(i) for i = 1:N]
end

# --------------------------------------------------------------------------------------------------
# Generatres an N_steps × N_T array of temperatures where N_T = length(temps) and each column k is a series
# of temperature steps that starts with T_start and ends in temps[k] and this series is geometrically
# distributed. Increase exponent to make the temperature steps more dense around lower temperatures.
function genGeometricTemperatureSteps(T_start::R, temps::Array{R, 1}, N_steps::Int; exponent=1) where R<:Real
    N_T = length(temps)
    T_mt = Array{R, 2}(undef, N_steps, N_T)
    for k = 1:N_T
        T_mt[:, k] = geometricSeries(T_start, temps[k], N_steps; exponent=exponent)
    end
    T_mt
end
function genGeometricTemperatureSteps(init_temps::Array{R, 1}, temps::Array{R, 1}, N_steps::Int; exponent=1) where R<:Real
    N_T = length(temps)
    T_mt = Array{R, 2}(undef, N_steps, N_T)
    for k = 1:N_T
        T_mt[:, k] = geometricSeries(init_temps[k], temps[k], N_steps; exponent=exponent)
    end
    T_mt
end


# --------------------------------------------------------------------------------------------------
# Returns a matrix of 2D momentum vectors in the 1BZ
function getMomMatrix(L::Int64)
    return [[2π/L*(x-1-L/2), 2π/L*(L/2-y)] for y=1:L, x=1:L]
end

# --------------------------------------------------------------------------------------------------
# Translates the 2d matrix V in each dimension by an index length tr_i and tr_j. Elements that fall
# outside the length of the array wraps around.
function translate2DMat(V::Array{T, 2}, tr_i::Int, tr_j::Int) where {T}
    Li = size(V, 1); Lj = size(V, 2)
    V_new = Array{T, 2}(undef, Li, Lj)
    for i = 1:Li, j = 1:Lj
        V_new[Int(mod(i + tr_i - 1, Li)+1), Int(mod(j + tr_j - 1, Lj)+1)] = V[i,j]
    end
    V_new
end

# --------------------------------------------------------------------------------------------------
# Given two matrices where the second is assumed to contain the errors of the matrix elements
# of the first, we calculate the relative error and returns the one with the largest relative error
# as well as its matrix indices.
function maxRelErr(avg_A::Array{T,2}, err_A::Array{T,2}) where T <: Real
    v = 1; h = 1
    val = err_A[1]/avg_A[1]
    for v_pos = 1:size(avg_A,1), h_pos = 1:size(avg_A,2)
        rel = err_A[v_pos,h_pos]/avg_A[v_pos,h_pos]
        if rel > val
            val = rel
            v = v_pos; h = h_pos
        end
    end
    return v, h, val
end


# --------------------------------------------------------------------------------------------------
# Takes in two numbers where the second is assumed to be the error. Rounds this to one significant
# digit and rounds the first value to the appropriate decimal place.
function scientificRounding(val::T, err::T; extra_digits::Int64 = 0) where T<:Real
    if err < 0.0
        println("Warning: error less than zero inserted: $(err)")
        err = abs(err)
    end
    # First we find the number of decimals needed for the first significant digit of the error
    # Round to first significant digit
    err_temp = round(err; sigdigits=1)#signif(err, 1)
    # Find the first significant digit
    st = "$(err_temp)"
    if st == "NaN"
        # Infinite error means we don't know what the value is.
        println("Warning: error was NaN")
        return signif(val, 1+extra_digits), err
    end
    if occursin(r"[e]", st)#ismatch(r"[e]", st)
        st = st[1]
    else
        st = reverse("$(err_temp)")
        st = match(r"[1-9]", st).match
    end
    sig = parse(Int64, st)
    # Now divide by this integer to get a number on the form 10^(-d) and extract d through log10
    digi_num = floor(Int64, -log10(err_temp/sig))
    # Finally use the number of digits to round the value (note that this number could be negative)
    val = round(val; digits=digi_num+extra_digits)
    return val, round(err; sigdigits=1+extra_digits)#signif(err, 1+extra_digits)
end

# --------------------------------------------------------------------------------------------------
# A series of matrices produced by Monte Carlo measurements, calculate the average and error
# for each element in the matrix.
function avgErr(A::Array{Array{Float64, 2},1})
    avg_A = mean(A)
    sm_A = zeros(avg_A)
    M = length(A)
    N = length(avg_A)
    for m = 1:M
        for i = 1:N
            sm_A[i] += A[m][i]^2
        end
    end
    sm_A = sm_A./M
    τ_matrix = [M/effective_sample_size([A[m][i] for m=1:M]) for i=1:N]
    τ_matrix = reshape(τ_matrix, size(avg_A)...)
    
    err_A = sqrt.((1+2 .*τ_matrix).*abs.(sm_A - avg_A.^2)./(M-1))
    
    return avg_A, err_A
end

# --------------------------------------------------------------------------------------------------
# A series of measurements produced by MC Monte Carlo calculations, calculate the average and error
# of the array.
function avgErr(A::Array{T,1}) where T<:Real
    avg_A = mean(A)
    sm_A = 0.0
    M = length(A)
    M > 1 || throw(error("ERROR: List of size >1 required"))
    for m = 1:M
        sm_A += A[m]^2
    end
    sm_A = sm_A/M
    τ = M/effective_sample_size(A)
    if τ<0
        println("Warning: Calculated correlation time is < 0. That is weird.")
        τ = abs(τ)
    end
    err_A = sqrt((1+2*τ)*abs(sm_A - avg_A^2)/(M-1))
    
    return avg_A, err_A
end

# --------------------------------------------------------------------------------------------------
# Produces a string suitable for showing the matrix element with maximum relative error.
function maxRelErrString(avg_A::Array{T,2}, err_A::Array{T,2}) where T <: Real
    v, h, rel_err = maxRelErr(avg_A, err_A)
    avg, err = scientificRounding(avg_A[v,h], err_A[v,h])
    return " rel_err = $(signif(rel_err,1))\t@ [$(v), $(h)] for\t $(avg)±$(err)"
end

# --------------------------------------------------------------------------------------------------
# Given a matrix with vertical dimension L_v, and the single number for an element in the matrix i,
# calculate the corresponding matrix indices v and h, by using the formula i = (h-1)*L_v+v
function matrixIndices(i::Int64, L_v::Int64)
    if i%L_v == 0
        v = L_v
        h = floor(Int64, i/L_v)
    else
        h = floor(Int64, i/L_v)+1
        v = i-(h-1)*L_v
    end
    return v,h
end

# --------------------------------------------------------------------------------------------------
# Find the n matrix elements with the highest value and return their matrix indices [v,h] in a list.
function findMaximaIndices(A::Array{T, 2}; n::Int64=1) where T <: Real
    L_v = size(A,1)
    maxes = [[A[i], i] for i = 1:length(A)]
    sort!(maxes, rev=true, by = x -> x[1])
    return [[matrixIndices(Int(maxes[i][2]), L_v)...] for i = 1:n]
end

# --------------------------------------------------------------------------------------------------
# Returns a sorted array of divisors of M except for 1 and M.
function divisors(M::Int; sorting=true)
    # First we calculate all the divisors of M except for 1 and M
    # and save them in the divisors array.
    pr = primes(2,Int(floor(sqrt(M)))+2)
    L = size(pr, 1)
    divisors = Array{Int}(undef, 0)
    m = M
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
    if sorting
        return sort(unique(divisors), rev=true)
    else
        return unique(divisors)
    end
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
    divs = divisors(M)
    for div in divs
        if div <= N
            return div
        end
    end
    
end

function mkcd(dir_name::AbstractString; visible=false)
    if ispath(dir_name)
        if visible
            println("Directory $(dir_name) already exists.. entering.")
        end
        cd(dir_name)
    else
        mkdir(dir_name)
        cd(dir_name)
    end
end

function meanAndErr(A::Vector{R}) where R<:Real
    mn = mean(A)
    er = std(A)/sqrt(length(A))
    (mn, er)
end

# Removing mid-point of L×L matrix
function removeMiddle(A::Array{T, 2}) where {T}
    L = size(A, 1)
    L_mid = floor(Int64, L/2)+1
    val = A[1,1]
    A′ = copy(A)
    A′[L_mid-1, L_mid] = val
    A′[L_mid-1, L_mid-1] = val
    A′[L_mid-1, L_mid+1] = val
    A′[L_mid-2, L_mid] = val
    A′[L_mid, L_mid] = val
    #A′[L_mid, L_mid+1] = val
    A′
end


