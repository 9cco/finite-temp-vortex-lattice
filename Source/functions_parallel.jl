using MCMCDiagnostics
using Base.Test

####################################################################################################
#                             Functions used for parallell measurements
#
####################################################################################################

# --------------------------------------------------------------------------------------------------
function nMCS(ψ::State, sim::Controls, n::Int64)
    for i=1:n
        mcSweep!(ψ, sim)
    end
    return ψ
end

# --------------------------------------------------------------------------------------------------
# Perform n MCsweeps and return the final state and an array with energies after each sweep
function nMCSEnergy(ψ::State, sim::Controls, n::Int64, E₀::Float64)
	E_array = Array{Float64,1}(n)
	E_array[1] = E₀ + mcSweepEn!(ψ,sim)
    for i=2:n
        E_array[i] = E_array[i-1] + mcSweepEn!(ψ,sim)
        #Have checked that this energy matches E(ψ) for each step
    end
    return ψ, E_array
end

#---------------------------------------------------------------------------------------------------
# Perform n MCsweeps and return the final state and an array with energies after each sweep
# and updates simulation constants dynamically
function nMCSEnergyDynamic(ψ::State, sim::Controls, n::Int64, adjust_int::Int64, E₀::Float64)
	Es = Array{Float64,1}(n)
    adjustment_mcs = 0

	Es[1] = E₀ + mcSweepEn!(ψ,sim)
    
    for i=2:n
        Es[i] = Es[i-1] + mcSweepEn!(ψ,sim)
        
        #Every adjust_int steps, update simulation constants if needed
        if i % adjust_int == 0
            #Testing
            (ar, mcs) = adjustSimConstants!(sim, ψ)
            adjustment_mcs += mcs
            
            #Update energy (since we do MCS in the update steps)
            Es[i] = E(ψ)
        end
    end
    return ψ, Es, sim, adjustment_mcs
end    

#--------------------------------------------------------------------------------------------------
# Cannot use adjustSimConstant in a paralell, process since it does not return ψ and sim, so have to
# make a dummy function to help send variables between processes.
function adjustSimConstantsPar(sim::Controls, ψ::State, M::Int64=40)
    accept_ratio, adjustment_mcs = adjustSimConstants!(sim, ψ, M)
    return ψ, sim, adjustment_mcs
end

function meanAndStd{T<:Real}(E₁::Array{T,1}, E₂::Array{T,1})
	dE = E₁ - E₂
	return mean(dE), std(dE)
end

#--------------------------------------------------------------------------------------------------
# Check if average from t_start to t_finish is the same in E_workers and E_ref where these
# arrays are 
function checkThermalization(E_ref::Array{Float64,1}, E_workers::Array{Float64,2}, n_workers::Int,
        t_start::Int64, t_end::Int64, STD_NUMBER::Float64=3.0)
    max_av = 0.0
    max_std = 0.0

	# Check that we have enough available processes
	nprocs()-1 >= n_workers || throw(error("ERROR: Not enough workers"))

	future_array = [Future() for i = 1:n_workers]
	av = Array{Float64,1}(n_workers)
	st = Array{Float64,1}(n_workers)

	for w = 1:n_workers
		future_array[w] = @spawn meanAndStd(E_ref[t_start:t_end], E_workers[w,t_start:t_end])
	end
	for w = 1:n_workers
		av[w], st[w] = fetch(future_array[w])
		av[w] = abs(av[w])
#		st[w] = st[w]/√(t_end-t_start+1) # Calculates standard error.
	end

	# Check if av <= std for all arrays and calculate the worker with the
	# maximum average
	thermalized = true
    therm_w = [true for i = 1:n_workers]
	for w=1:n_workers
		if av[w] > st[w]*STD_NUMBER
			thermalized = false
            therm_w[w] = false
		elseif av[w] > max_av
			max_av = av[w]
			max_std = st[w]
		end
	end

    if thermalized
		println("Thermalization successful between T = [$(t_start), $(t_end)]")
		for w = 1:n_workers
			println("Worker $(w): ΔE = $(av[w]) ± $(st[w])")
		end
        return true
    else
        for w=1:n_workers
            if therm_w[w] == false
                println("Problem for worker $(w): ΔE = $(av[w]) ± $(st[w])")
                break
            end
        end
        return false
    end
end

# --------------------------------------------------------------------------------------------------
# We assume that ψ has already reached equilibirum
function parallelMultiplyState(ψ::State, sim::Controls, t₀::Int64)
    # Find the number of necessary equilibrium states
    np = nprocs()-1 # One of the processes is the master process.
    # Copy state to list until we have a storage for all the necessary worker states.
    ψ_list = [copy(ψ) for i=1:np+1]
    future_list = [Future() for i=1:np]
    
    # Start np workers doing t₀ number of Monte-Carlo Sweeps on ψ
    for i=1:np
        future_list[i] = @spawn nMCS(ψ, sim, t₀)
    end
    
    # Use master process to create a final state
    ψ_list[np+1] = nMCS(ψ, sim, t₀)
    
    # Now save the worker states into the list
    for i=1:np
        ψ_list[i] = fetch(future_list[i])
    end
    
    return ψ_list
end

# --------------------------------------------------------------------------------------------------
# Makes M measurements of both the vortex lattice and the structure factor.
function sfvlaMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64, Δt::Int64)
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
    for m = 2:M
#        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
#        flush(STDOUT)
        
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
function parallelSFVLA!{T<:Real}(ks::Array{Array{T, 1}, 2}, 
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
on a $(L)×$(L) system.")
    
    # Start +1 workers
    for i = 1:nw
        futures[i] = @spawn sfvlaMeasure!(ks, ψ_list[i], sim, M_min+1, Δt)
    end
    # Start remaining workers
    for i = 1:np-nw-1
        futures[nw+i] = @spawn sfvlaMeasure!(ks, ψ_list[nw+i], sim, M_min, Δt)
    end
    # Make the master process work as well
    (S⁺, S⁻, V⁺, V⁻) = sfvlaMeasure!(ks, ψ_list[np], sim, M_min, Δt, "-v")
    
    println("Measurements done, collecting parallell results.")
    # Collect results
    for i = 1:np-1
        (new_S⁺, new_S⁻, new_V⁺, new_V⁻) = fetch(futures[i])
        S⁺ = vcat(S⁺, new_S⁺)
        S⁻ = vcat(S⁻, new_S⁻)
        V⁺ = vcat(V⁺, new_V⁺)
        V⁻ = vcat(V⁻, new_V⁻)
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
    
    return (av_V⁺, err_V⁺, V⁺, av_V⁻, err_V⁻, V⁻, av_S⁺, err_S⁺, S⁺, av_S⁻, err_S⁻, S⁻)
end

# Find first occurrence of when the values in the arrays cross eachother in the interval
# between iₛ and iₑ
function firstZero{T<:Real}(E_high::Array{T, 1}, E_low::Array{T,1}, iₛ::Int64, iₑ::Int64)
	for i = iₛ:iₑ
		if E_high[i] <= E_low[i]
			return i
		end
	end
	return -1
end

#--------------------------------------------------------------------------------------------------
# Thermalises several states simultaneously on multiple cores. One high T state is thermalised
# along with (#available cores - 1) low T states, where energies are compared to ensure proper
# convergence.
function parallelThermalization!(ψ_ref::State, ψ_w::Array{State,1}, c::SystConstants,
        sim::Controls, T::Int64=1000, ex::Float64=1.8, STD_NUMBER::Float64=0.5)
    CUTOFF_MAX::Int64=4000000           #Max number of MCS before the function terminates
    ADJUST_INTERVAL=400                 #Number of MCS between each sim_const adjustment while finding dE<0
    AVERAGING_INT_FRAC=1/4                 #Similiar for 2nd loop, also the interval that is averaged over
    AVERAGING_INT_EX = 1.4

    #Initialisation
    NWS=length(ψ_w)                     #Number of states in ψ_list
	E_w = Array{Float64, 2}(NWS, CUTOFF_MAX)        #Matrix to store energies of worker states
	for w = 1:NWS
		E_w[w, 1] = E(ψ_w[w])
	end
	E_ref = Array{Float64, 1}(CUTOFF_MAX)           #Array to store energies of reference state
	E_ref[1] = E(ψ_ref)

	# Checking wheather reference or worker has highest energy.
	ref_highest = [true for i = 1:NWS]
	for w = 1:NWS
		if E_w[w,1]>E_ref[1]
			ref_highest[w] = false
		end
	end

	T += 1								# T is as a parameter the initial number of steps taken to find the thermalization and
										# additionally, in the loops it is the time-index a loop should reach.

    sim_ref = copy(sim)                 #Simulation constants for reference state(s)
    sim_w = [copy(sim) for i=1:NWS]       #Simulation constants for the workers
    adjustment_mcs = 0                  #Number of MCsweeps done during dynamic updates sysconst
    mcs_ref = 0                         # Contains number of MCS done in adjustments of ref state
    mcs_w = zeros(Int64, NWS)                  # Contains number of MCS done in adjustments of worker states.
    E_check_workers = zeros(Int64, NWS) #Check the requirement dE <= 0 for each process
    tₛ=2                                # First free time-index after last update of energy array.
    thermalized_init = false
    N = (ψ_ref.consts.L)^2              #Number of lattice sites

    #Find the number of available workers
    np = nprocs()-1
    println("Number of parallel workers: $(np)")
    #Check if we have more low T states than workers
    if NWS > np
        error("Not enough available workers, $(np) workers and $(NWS) states")
    end
    
    # Make a future list to use for parallel thermalisation
    ψ_future_list = [Future() for i=1:np]
    
    # Intial thermalisation until the high T state has lower energy than low T at some point
    while !thermalized_init && T < CUTOFF_MAX 
        #Run T MCsweeps for all the workers with dynamic sysconst updates
        for w=1:NWS
			ψ_future_list[w] = @spawn nMCSEnergyDynamic(ψ_w[w], sim_w[w], T-tₛ+1, ADJUST_INTERVAL, E_w[w,tₛ-1])
        end
        #Similar for the high T in master
		ψ_ref, E_ref[tₛ:T], sim_ref, mcs_ref = nMCSEnergyDynamic(ψ_ref, sim_ref, T-tₛ+1, ADJUST_INTERVAL, E_ref[tₛ-1])
        
        #Gather results from workers
        for w=1:NWS
            ψ_w[w], E_w[w,tₛ:T], sim_w[w], mcs_w[w] = fetch(ψ_future_list[w])
        end
        adjustment_mcs += max(maximum(mcs_w),mcs_ref)
        
        # Check if for all processes dE < 0, and get largest step that this happens for.
		# Takes into account whether it is the worker or reference that should have highest energy.
		for w = 1:NWS
			if ref_highest[w]
				ψ_future_list[w] = @spawn firstZero(E_ref, E_w[w, :], tₛ, T)
			else
				ψ_future_list[w] = @spawn firstZero(E_w[w,:], E_ref, tₛ, T)
			end
		end
		for w = 1:NWS
			i = fetch(ψ_future_list[w])
			if i != -1
                E_check_workers[w] = i
                println("Worker $(w) initially thermalised after $(i) steps")
                flush(STDOUT)
			end
		end

        if minimum(E_check_workers) == 0 #This is always true if we found no dE < 0 in all workers
            tₛ = T
            T = min(Int(ceil(T*ex)), CUTOFF_MAX)
			println("Increasing simulation time such that tₛ = $(tₛ) and T = $(T)")
			flush(STDOUT)
            if T == CUTOFF_MAX
                println("Failed to find a point where ΔE <= 0")
                return (-1, tₛ, T,  E_ref[1:tₛ], E_w[:,1:tₛ], ψ_ref, ψ_w, sim_ref, sim_w)
            end
        else
            tₛ = maximum(E_check_workers)
            println("All workers initially thermalized after $(T) steps")
			flush(STDOUT)
            thermalized_init = true
        end
    end

    tₛ = T+1
	averaging_int = ceil(Int64, (T+adjustment_mcs)*AVERAGING_INT_FRAC)
    T = tₛ + averaging_int - 1
    #Set this to A_I so that the simulation constants are only updated by energyDynamic at the END of each 
    #iteration in the while loop.
    #Increase simulation time by ADJUST_INT_THERM and try to find an interval with small average dE
    while T<CUTOFF_MAX
		println("Checking average ∈ [$tₛ, $T]") 

        #Start by updating the simulation constants
        for w=1:NWS
            ψ_future_list[w] = @spawn(adjustSimConstantsPar(sim_w[w], ψ_w[w]))
        end
        ψ_ref, sim_ref, mcs_ref = adjustSimConstantsPar(sim_ref, ψ_ref)
        for w=1:NWS
            ψ_w[w], sim_w[w], mcs_w[w] = fetch(ψ_future_list[w]) 
        end
        adjustment_mcs += max(maximum(mcs_w),mcs_ref)

		# Correct energy for built-up floating point errors
		E_ref[tₛ-1] = E(ψ_ref)
		for w=1:NWS
			E_w[w,tₛ-1] = E(ψ_w[w])
		end
        
        # Do MCS to find average
        for w=1:NWS
			ψ_future_list[w] = @spawn nMCSEnergy(ψ_w[w], sim_w[w], averaging_int, E_w[w,tₛ-1])
        end
		ψ_ref, E_ref[tₛ:T] = nMCSEnergy(ψ_ref, sim_ref, averaging_int, E_ref[tₛ-1])
        for w=1:NWS
            ψ_w[w], E_w[w,tₛ:T] = fetch(ψ_future_list[w])
        end

        if checkThermalization(E_ref./N, E_w./N, NWS, tₛ, T,STD_NUMBER)
			println("Final thermalization time: $(T+adjustment_mcs)")
            return(T+adjustment_mcs, tₛ, T, E_ref[1:T], E_w[:,1:T], ψ_ref, ψ_w, sim_ref, sim_w)
        end
        
        # If we didn't find thermalization in this interval, we extend the interval by a fraction and go
        # to this new interval.
        averaging_int = ceil(Int64, averaging_int*AVERAGING_INT_EX)
        tₛ=T + 1
        T = tₛ + averaging_int - 1
        println("Increasing simulation time to T = $(T)")
    end
    return -1, tₛ, T, E_ref, E_w, ψ_ref, ψ_w, sim_ref, sim_w
end
