###################################################################################################
#               Functions to measure the order parameters, along with vorticity etc
###################################################################################################
# Makes M measurements of both the vortex lattice and the structure factor.
function orderParamMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64, Δt::Int64)
    L = ψ.consts.L
    L_k = size(ks, 1)
    # Setup measurement storage
    S⁺ = [zeros(L_k,L_k) for i=1:M]
    S⁻ = [zeros(L_k,L_k) for i=1:M]
    V⁺ = [zeros(L,L) for i=1:M]
    V⁻ = [zeros(L,L) for i=1:M]
    
    Z₂ = [zeros(Complex, L_k, L_k) for i=1:M]
    Z₂_local_temp = Array{Complex}(L,L)
    GS = Array{Float64}(M)
    k_GS = [2*π/L, 0.0] #For which k-value we wish to measure Gausestiffness
    mid_v = Int(ceil(L/2))
    mid_h = Int(ceil(L/2))
    
    # The first measurement is of the initial ψ
    (V⁺[1], V⁻[1]) = vortexSnapshot(ψ)
    
    #Measure order parameters    
    GS[1] = gaugeStiffnessMeasure(k_GS, ψ)
    Z₂MeasureLocal!(ψ, Z₂_local_temp)

    # Then we use this to measure the structure factor
    for x=1:L_k, y=1:L_k
        (S⁺[1][y,x], S⁻[1][y,x]) = structureFunction(ks[y,x], ψ, V⁺[1], V⁻[1])
        Z₂[1][y,x] = Z₂StructureFunction(ks[y,x], ψ, Z₂_local_temp, mid_h, mid_v, L)
    end
    
    # For each measurement
    for m = 2:M
#        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
#        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Find n_z(r) of the lattice.
        (V⁺[m], V⁻[m]) = vortexSnapshot(ψ)
        #Measure order parameters    
        GS[m] = gaugeStiffnessMeasure(k_GS, ψ)
        Z₂MeasureLocal!(ψ, Z₂_local_temp)


        # Find structure factor. 
        for x=1:L_k, y=1:L_k
            (S⁺[m][y,x], S⁻[m][y,x]) = structureFunction(ks[y,x], ψ, V⁺[m], V⁻[m])
            Z₂[m][y,x] = Z₂StructureFunction(ks[y,x], ψ, Z₂_local_temp, mid_h, mid_v, L)
        end
    end
    
    # After the loop, we should have filled up the M measurements for the matrices.
    return (S⁺, S⁻, V⁺, V⁻, Z₂, GS)
end

# Same as above, but now prints progress to STDOUT
function sfvlaMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64,
        Δt::Int64, option::AbstractString)
	PROG_NUM = 10		# Adjusts the number of times progress is reported while measuring.
    L = ψ.consts.L
    L_k = size(ks, 1)
    # Setup measurement storage
    S⁺ = [zeros(L_k,L_k) for i=1:M]
    S⁻ = [zeros(L_k,L_k) for i=1:M]
    V⁺ = [zeros(L,L) for i=1:M]
    V⁻ = [zeros(L,L) for i=1:M]
    
    # The first measurement is of the initial ψ
    (V⁺[1], V⁻[1]) = vortexSnapshot(ψ)
    
    # Then we use this to measure the structure factor
    for x=1:L_k, y=1:L_k
        (S⁺[1][y,x], S⁻[1][y,x]) = structureFunction(ks[y,x], ψ, V⁺[1], V⁻[1])
    end
    
    # For each measurement
	prog_int = floor(Int64, M/PROG_NUM)
    for m = 2:M
		if m % prog_int == 0
			println("Measurement progress: $(Int(round(m/M*100,0)))%")
        	flush(STDOUT)
		end
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Find n_z(r) of the lattice.
        (V⁺[m], V⁻[m]) = vortexSnapshot(ψ)
        
        # Find structure factor. 
        for x=1:L_k, y=1:L_k
            (S⁺[m][y,x], S⁻[m][y,x]) = structureFunction(ks[y,x], ψ, V⁺[m], V⁻[m])
        end
    end
    
    # After the loop, we should have filled up the M measurements for the matrices.
    return (S⁺, S⁻, V⁺, V⁻)
end

# --------------------------------------------------------------------------------------------------
# Given np+1 uncorrelated states in ψ_list we use these to make M measurements of the vortex lattice by splitting
# the M measurements on the np workers as well as the master process. In this version we continuously discard
# the measurements and only save the averages and second moments.
function parallelMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, 
        ψ_list::Array{State,1}, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ_list[1].consts
    L = syst.L
    
    # Checking that the k matrix has equal dimensions
    L_k = size(ks, 1)
    size(ks, 2) == L_k || throw(DomainError())
    
    s_norm_inv = 1/(L^2*syst.f*two_pi)^2
    v_norm_inv = 1/two_pi
    
    # Splitting the problem into np sub-problems.
    np = nprocs()
    # Minimum amount of work pr. process
    M_min = Int(floor(M/np))
    # Number of workers doing +1 extra work
    nw = M%np
    
    # Make sure that we have enough states
    length(ψ_list) >= np || throw(error("ERROR: Not enough states in list"))
    
    # Setup worker futures
    futures = [Future() for i=1:(np-1)]
    
    println("Starting $(M) measurements on $(np) processes doing max $(M_min + Int(ceil(nw/np))) measurements each
on a $(L)×$(L) system, corresponding to $((M_min+ceil(Int64, nw/np))*Δt) MCS pr. process")
    
    # Start +1 workers
    for i = 1:nw
        futures[i] = @spawn orderParamMeasure!(ks, ψ_list[i], sim, M_min+1, Δt)
    end
    # Start remaining workers
    for i = 1:np-nw-1
        futures[nw+i] = @spawn orderParamMeasure!(ks, ψ_list[nw+i], sim, M_min, Δt)
    end
    # Make the master process work as well
    (S⁺, S⁻, V⁺, V⁻, Z₂, GS) = orderParamMeasure!(ks, ψ_list[np], sim, M_min, Δt)
    
    println("Measurements done, collecting parallell results.")
    # Collect results
    for i = 1:np-1
        (new_S⁺, new_S⁻, new_V⁺, new_V⁻, new_Z₂, new_GS) = fetch(futures[i])
        S⁺ = vcat(S⁺, new_S⁺)
        S⁻ = vcat(S⁻, new_S⁻)
        V⁺ = vcat(V⁺, new_V⁺)
        V⁻ = vcat(V⁻, new_V⁻)
        Z₂ = vcat(Z₂, new_Z₂)
        GS = vcat(GS, new_GS)
    end
    @test length(S⁺) == M
    
    println("Parallell measurements done. Processing.")
    
    # Normalizing results
    V⁺ = v_norm_inv.*V⁺
    V⁻ = v_norm_inv.*V⁻
    S⁺ = s_norm_inv.*S⁺
    S⁻ = s_norm_inv.*S⁻
    
    # Calculate averages and second moments
    av_S⁺ = mean(S⁺)
    av_S⁻ = mean(S⁻)
    av_V⁺ = mean(V⁺)
    av_V⁻ = mean(V⁻)
    av_Z₂ = mean(Z₂)
    av_GS = mean(GS)
    sm_S⁺ = zeros(L_k, L_k)
    sm_S⁻ = zeros(L_k, L_k)
    sm_V⁺ = zeros(L,L)
    sm_V⁻ = zeros(L,L)
    for m = 1:M
        sm_S⁺ += S⁺[m].^2
        sm_S⁻ += S⁻[m].^2
        sm_V⁺ += V⁺[m].^2
        sm_V⁻ += V⁻[m].^2
    end
    sm_S⁺ = sm_S⁺./M
    sm_S⁻ = sm_S⁻./M;
    sm_V⁺ = sm_V⁺./M
    sm_V⁻ = sm_V⁻./M
    
    # Error calculation of average vorticity
    τ_V⁺ = [1/ess_factor_estimate([V⁺[m][y,x] for m=1:M])[1] for y=1:L, x=1:L]
    τ_V⁻ = [1/ess_factor_estimate([V⁻[m][y,x] for m=1:M])[1] for y=1:L, x=1:L]
    
    err_V⁺ = (1+2.*τ_V⁺).*(sm_V⁺ - av_V⁺.^2)./(M-1)
    err_V⁻ = (1+2.*τ_V⁻).*(sm_V⁻ - av_V⁻.^2)./(M-1)
    
    # Error calculation of structure factor.
    τ_S⁺ = [1/ess_factor_estimate([S⁺[m][y,x] for m=1:M])[1] for y=1:L_k, x=1:L_k]
    τ_S⁻ = [1/ess_factor_estimate([S⁻[m][y,x] for m=1:M])[1] for y=1:L_k, x=1:L_k]
    
    err_S⁺ = (1+2.*τ_S⁺).*(sm_S⁺ - av_S⁺.^2)./(M-1)
    err_S⁻ = (1+2.*τ_S⁻).*(sm_S⁻ - av_S⁻.^2)./(M-1)

	# Sum of all vorticities
	println("\nSum of vorticity of random snapshot:\nV⁺: \t$(sum(V⁺[rand(1:M)]))
V⁻: \t$(sum(V⁻[rand(1:M)]))")
    
    
    # Finding max values over the matrices.
    max_S⁺= maximum(av_S⁺)
    max_S⁻ = maximum(av_S⁻)
    max_err_S⁺ = maximum(err_S⁺)
    max_err_S⁻ = maximum(err_S⁻)
    max_τ_S⁺ = maximum(τ_S⁺)
    max_τ_S⁻ = maximum(τ_S⁻)
    println("\nMax (S⁺, S⁻)\n($(max_S⁺), $(max_S⁻))")
    println("Max δ(S⁺, S⁻)\n($(max_err_S⁺), $(max_err_S⁻))")
    println("Max correlation time\n($(max_τ_S⁺), $(max_τ_S⁻))")
    
    return (av_V⁺, err_V⁺, V⁺, av_V⁻, err_V⁻, V⁻, av_S⁺, err_S⁺, S⁺, av_S⁻, err_S⁻, S⁻, av_Z₂, av_GS)
end



