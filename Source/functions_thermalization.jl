####################################################################################################
#                             Thermalization functions
#
####################################################################################################


# ---------------------------------------------------------------------------------------------------
# Version of find Equilibrium that doesn't waste MCS when adjusting sim Constants and also calculates energy based
# on return value of mcSweepEn! as well as gets rid of the use of dynamic arrays.
# ---------------------------------------------------------------------------------------------------
# Given values for the physical constants of the system as well as the system size, we find the number of MC-sweeps it
# takes until the internal energy of the system reaches a more or less constant value.
function findEquilibrium(c::SystConstants, sim₁::Controls=Controls(π/3, 0.4, 3.0), 
        T::Int64=1000, ex::Float64=1.5, di::Int64=8)
    CUTOFF_MAX::Int64 = 8000000
    ADJUST_INTERVAL::Int64 = 2000
    STD_NUMBER::Int64 = 1
    println("Finding Equilibrium of\n$(c)\n$(sim₁)")
    
    ψ₂ = State(2, c)
    sim₂ = copy(sim₁)
    ψ₁ = State(1, c)
    dE = zeros(CUTOFF_MAX)
    E₁ = zeros(CUTOFF_MAX)
    E₂ = zeros(CUTOFF_MAX)
    adjustment_mcs = 0
    
    # Check that the un-correllated state has higher energy than the correlated
    if E(ψ₂) <= E(ψ₁)
        error("Correlated state has higher energy than un-correlated")
    end
    
    tₛ = 0    # The wanted t₀ does not exist at or before this position.
    t₀ = T

	println("Performing initial $(T) MCS")
	flush(STDOUT)
    
    E₁[1] = E(ψ₁)
    E₂[1] = E(ψ₂)
    dE[1] = E₂[1] - E₁[1]
    for i = 2:T
        E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
        E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
        dE[i] = E₂[i] - E₁[i]
    end
    
    # Adjust simulation constants as needed
    (ar, mcs1) = adjustSimConstants!(sim₁, ψ₁)
    (ar, mcs2) = adjustSimConstants!(sim₂, ψ₂)
    adjustment_mcs += max(mcs1, mcs2)   # Adds the max number of monte-carlo sweeps done for the two
    # adjustments to the number of adjustments used for finding the equilibrium time.
    
    E₁[T] = E(ψ₁)
    E₂[T] = E(ψ₂)
    dE[T] = E₂[T] - E₁[T]
    
    while tₛ < CUTOFF_MAX
        # Find the first occurence of dE <= 0 if it exists
        println("Searching for ΔE <= 0..")
		flush(STDOUT)
        t₀ = T
        for i = (tₛ+1):T
            if dE[i] <= 0
                t₀ = i
                break
            end
        end
        
        while T <= t₀ && T < CUTOFF_MAX
            # If we couldn't find a t₀ in dE we have to try and increase simulation time
            tₛ = T
            T = min(Int(ceil(T*ex)), CUTOFF_MAX)
            for i = (tₛ+1):T
                E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
                E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
                dE[i] = E₂[i] - E₁[i]
                
                # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
                if i % ADJUST_INTERVAL == 0
                    ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                    ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                    adjustment_mcs += max(mcs1, mcs2)
                    
                    # Then we also go over an extra time to get energy correct
                    E₁[i] = E(ψ₁)
                    E₂[i] = E(ψ₂)
                    dE[i] = E₂[i] - E₁[i]
                end
            end
            
            # Then we again see if we can find the first occurrence of dE <= 0 after tₛ
            t₀ = T
            for i = tₛ:T
                if dE[i] <= 0
                    t₀ = i
                    break
                end
            end
            
            if t₀ == T == CUTOFF_MAX # We have not found any dE < 0 and we have reached the max number of sweeps
                println("Failed to find a point where ΔE <= 0")
                return (-1, E₁, E₂, dE, ψ₁, ψ₂, sim₁, sim₂)
            end
            
            # When the loop ends we should have the situation that T > t₀ where t₀ is the first occurrence
            # in dE where dE[t₀] <= 0
        end
        println("ΔE <= 0 found at t₀ = $(t₀)!\nChecking if average is close to 0..")
		println("$(Int(round(T/CUTOFF_MAX*100,0)))% of max") # Debug
		flush(STDOUT)
        
        # Now we make sure that T is large enough such that [1,T] includes an interval [t₀, t₀+t₀/div]
        # so that an average can be performed
        t_end = min(t₀ + Int(ceil(t₀/di)), CUTOFF_MAX)
        while T < t_end
            T += 1
            E₁[T] = E₁[T-1] + mcSweepEn!(ψ₁, sim₁)
            E₂[T] = E₂[T-1] + mcSweepEn!(ψ₂, sim₂)
            dE[T] = E₂[T] - E₁[T]

            # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
            if T % ADJUST_INTERVAL == 0
                ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                adjustment_mcs += max(mcs1, mcs2)
                
                # Then we also go over an extra time to get energy correct
                E₁[T] = E(ψ₁)
                E₂[T] = E(ψ₂)
                dE[T] = E₂[T] - E₁[T]
            end
        end
        
        # Now we calculate the average and standard deviation of dE over [t₀, T] and check if the
        # average is within STD_NUMBER of standard deviations of 0 at which point we declare the equilibrium
        # time to have been found.
        int = dE[t₀:T]
        av = mean(int)
        st = std(int)
        if abs(av) <= STD_NUMBER*st
            println("Equilibrium found at time $(T+adjustment_mcs)
over the interval [$(t₀), $(T)]
s.t. <ΔE> = $(round(av,2)) ± $(round(st/sqrt(size(int,1)), 1))
std(ΔE) = $(round(st, 1))")
            return (T+adjustment_mcs, E₁[1:T], E₂[1:T], dE[1:T], ψ₁, ψ₂, sim₁, sim₂)
        end
        
        println("Average was not close to 0. Increasing interval.")
        
        # If we didn't find an interval that had an average close to 0 we assume this interval is ahead of us
        # and start again with an increased T, setting the starting point tₛ to the end of the interval.
        tₛ = T
        T = min(Int(ceil(T*ex)), CUTOFF_MAX)
        
        # Simulating new MCS
        for i = (tₛ+1):T
            E₁[i] = E₁[i-1] + mcSweepEn!(ψ₁, sim₁)
            E₂[i] = E₂[i-1] + mcSweepEn!(ψ₂, sim₂)
            dE[i] = E₂[i] - E₁[i]

            # After ADJUST_INTERVAL # MCS we see if adjusting the simulations constants is neccessary.
            if i % ADJUST_INTERVAL == 0
                ar, mcs1 = adjustSimConstants!(sim₁, ψ₁)
                ar, mcs2 = adjustSimConstants!(sim₂, ψ₂)
                adjustment_mcs += max(mcs1, mcs2)
                
                # Then we also go over an extra time to get energy correct
                E₁[i] = E(ψ₁)
                E₂[i] = E(ψ₂)
                dE[i] = E₂[i] - E₁[i]
            end
        end
    end
    return (-1, E₁, E₂, dE, ψ₁, ψ₂, sim₁, sim₂)
end
# IDEA: We could have sim be controlled as a Metropolis Update where the internal energy equivalent could be the number of
# accepts pr proposal over mcSweep and how close that is to 1/2.

#--------------------------------------------------------------------------------------------------
function meanAndStd{T<:Real, I<:Int}(E_ref::Array{T,1}, E_w::Array{T,2}, w::I, t_start::I, t_end::I)
    N = t_end-t_start+1
    dE = Array{T,1}(N)
    for i=1:N
        dE[i] = E_ref[t_start+i-1] - E_w[w,t_start+i-1]
    end
	return mean(dE), std(dE)
end

#--------------------------------------------------------------------------------------------------
# Check if average from t_start to t_finish is the same in E_workers and E_ref where these
# arrays are 
function checkThermalization(E_ref::Array{Float64,1}, E_workers::Array{Float64,2},
        t_start::Int64, t_end::Int64, STD_NUMBER::Float64=3.0)
    max_av = 0.0
    max_std = 0.0
    n_workers = size(E_workers,1)

	future_array = [Future() for i = 1:n_workers]
	av = Array{Float64,1}(n_workers)
	st = Array{Float64,1}(n_workers)

	for w = 1:n_workers
		#future_array[w] = @spawn meanAndStd(E_ref[t_start:t_end], E_workers[w,t_start:t_end])
		future_array[w] = @spawn meanAndStd(E_ref, E_workers, w, t_start, t_end)
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

#--------------------------------------------------------------------------------------------------
# Find first occurrence of when the values in the arrays cross eachother in the interval
# between iₛ and iₑ
function firstZero{T<:Real}(E_ref::Array{T, 1}, E_w::Array{T, 2}, w::Int64, iₛ::Int64, iₑ::Int64, ref_highest::Bool)
    if ref_highest
        for i = iₛ:iₑ
            if E_ref[i] <= E_w[w,i]
                return i
            end
        end
        return -1
    else
        for i = iₛ:iₑ
            if E_w[w,i] <= E_ref[i]
                return i
            end
        end
        return -1
    end
end

#--------------------------------------------------------------------------------------------------
# Thermalises several states simultaneously on multiple cores. One high T state is thermalised
# along with (#available cores - 1) low T states, where energies are compared to ensure proper
# convergence.
function parallelThermalization!(ψ_ref::State, ψ_w::Array{State,1}, c::SystConstants,
        sim::Controls, T::Int64=1000, ex::Float64=1.8, STD_NUMBER::Float64=0.5)
    CUTOFF_MAX::Int64=2000000           #Max number of MCS before the function terminates
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

	# Checking whether reference or worker has highest energy.
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
        println("WARNING: More requested states ($(NWS+1)), than currently available extra processes ($(np)).
Thermalization will be ×$(floor(Int64, (NWS+1)/(np+1))) as long.")
    end
    
    # Make a future list to use for parallel thermalisation
    ψ_future_list = [Future() for i=1:NWS]
    
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
            ψ_future_list[w] = @spawn firstZero(E_ref, E_w, w, tₛ, T, ref_highest[w])
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
        gc() # Trying to run garbage collection manually since Vilje seems to be running out of memory at long thermalizations.
    end

    tₛ = T+1
	averaging_int = ceil(Int64, (T+adjustment_mcs)*AVERAGING_INT_FRAC)
    T = tₛ + averaging_int - 1
    #Set this to A_I so that the simulation constants are only updated by energyDynamic at the END of each 
    #iteration in the while loop.
    #Increase simulation time by ADJUST_INT_THERM and try to find an interval with small average dE
    while T<CUTOFF_MAX
		println("Checking average ∈ [$tₛ, $T]")
        flush(STDOUT)

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

        if checkThermalization(E_ref, E_w, tₛ, T,STD_NUMBER)
			println("Final thermalization time: $(T+adjustment_mcs)")
            return(T+adjustment_mcs, tₛ, T, E_ref[1:T], E_w[:,1:T], ψ_ref, ψ_w, sim_ref, sim_w)
        end
        
        # If we didn't find thermalization in this interval, we extend the interval by a fraction and go
        # to this new interval.
        averaging_int = ceil(Int64, averaging_int*AVERAGING_INT_EX)
        tₛ = T + 1
        T = tₛ + averaging_int - 1
        println("Increasing simulation time to T = $(T)")
        gc() # Trying to run garbage collection manually since Vilje seems to be running out of memory at long thermalizations.
    end
    return -1, tₛ, CUTOFF_MAX, E_ref, E_w, ψ_ref, ψ_w, sim_ref, sim_w
end




