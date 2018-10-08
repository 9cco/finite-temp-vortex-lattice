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
#Perform n MCsweeps and return the final state and an array with energies after each sweep
function nMCS_energy(ψ::State, sim::Controls, n::Int64)
    E_array = zeros(n+1)
    E_array[1] = E(ψ)
    for i=1:n
        E_array[i+1] = E_array[i] + mcSweepEn!(ψ,sim)
        #Have checked that this energy matches E(ψ) for each step
    end
    return ψ, E_array
end
#---------------------------------------------------------------------------------------------------
#Same as nMCS_energy, but updates simulation constants dynamically
function nMCS_energyDynamic(ψ::State, sim::Controls, n::Int64, adjust_int::Int64)
    E_array = zeros(n+1)
    E_array[1] = E(ψ)
    adjustment_mcs = 0
    
    for i=1:n
        E_array[i+1] = E_array[i] + mcSweepEn!(ψ,sim)
        
        #Every adjust_int steps, update simulation constants if needed
        if i % adjust_int == 0
            #Testing
            (ar, mcs) = adjustSimConstants!(sim, ψ)
            adjustment_mcs += mcs
            
            #Update energy (since we do MCS in the update steps)
            E_array[i] = E(ψ)
        end
    end
    return ψ, E_array, sim, adjustment_mcs
end    

#--------------------------------------------------------------------------------------------------
#Cannot use adjustSimConstant in a paralell, process since it does not return ψ and sim, so have to
#make a dummy function to help send variables between processes.
function adjustSimConstantsPar(sim::Controls, ψ::State, M::Int64=40)
    accept_ratio, adjustment_mcs = adjustSimConstants!(sim, ψ, M)
    return ψ, sim, adjustment_mcs
end
#--------------------------------------------------------------------------------------------------
#Check if average from t_start to t_finnish is the same in E_workers and E_ref
function check_thermalisation(E_ref::Array{Float64,1}, E_workers::Array{Float64,2}, Nworkers::Int,
        t_start::Int64, t_end::Int64, STD_NUMBER::Float64=0.1)
    max_av = 0.0
    max_std = 0.0
    for i=1:Nworkers
        av = abs(mean(E_ref[t_start:t_end] - E_workers[i,t_start:t_end]))
        if av > max_av
            max_av = av
            max_std = std(E_ref[t_start:t_end] - E_workers[i,t_start:t_end])
        end
    end 
    if max_av <= STD_NUMBER*max_std
        return true
    else
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
		this_pr = Int(round(m/M*100,0))
		if this_pr % 10 == 0
        	println("Measurement progress: $(this_pr)%")
		end
        flush(STDOUT)
        
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
    length(ψ_list) >= np || throw(Error("ERROR: Not enough states in list"))
    
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

#--------------------------------------------------------------------------------------------------
#Thermalises several states simultaneously on multiple cores. One high T state is thermalised
# along with (#available cores - 1) low T states, where energies are compared to ensure proper
# convergence.
function parallelThermalisation!(ψ_ref::State, ψ_w::Array{State,1}, c::SystConstants,
        sim::Controls, T::Int64=1000, ex::Float64=1.5, STD_NUMBER::Float64=0.5)
    CUTOFF_MAX::Int64=800000            #Max number of MCS before the function terminates
    ADJUST_INTERVAL=400                 #Number of MCS between each sim_const update in 1st loop
    AVG_MCS_INTERVAL=4000               #Similiar for 2nd loop, also the interval that is averaged over

    #Initialisation
    NWS=length(ψ_w)                     #Number of states in ψ_list
    E_w = zeros(NWS, CUTOFF_MAX)        #Matrix to store energies of worker states
    E_ref = zeros(CUTOFF_MAX)           #Array to store energies of reference state
    sim_ref = copy(sim)                 #Simulation constants for reference state(s)
    sim_w = [copy(sim) for i=1:3]       #Simulation constants for the workers
    adjustment_mcs = 0                  #Number of MCsweeps done during dynamic updates sysconst
    mcs_ref = 0                         #Temp variable
    mcs_w = zeros(NWS)                  #Temp variables
    E_check_workers = zeros(Int64, NWS) #Check the requirement dE <= 0 for each process
    t₀=0                                #Keep track of the first while loop
    tₛ=0                                #Keep track of the second while loop
    thermalised_init = false
    N = (ψ_ref.consts.L)^2              #Number of lattice sites

    #Find the number of available workers
    np = nprocs()-1
    println("np = $(np)")
    #Check if we have more low T states than workers
    if NWS > np
        error("Not enough available workers, $(np) workers and $(NWS) states")
    end
    
    #Make a future list to use for paralell thermalisation
    ψ_future_list = [Future() for i=1:np]
    
    #Intial thermalisation until the high T state has lower energy than low T at some point
    while !thermalised_init && T < CUTOFF_MAX 
        #Run T MCsweeps for all the workers with dynamic sysconst updates
        for i=1:NWS
            ψ_future_list[i] = @spawn nMCS_energyDynamic(ψ_w[i], sim_w[i], T-t₀, ADJUST_INTERVAL)
        end
        #Similar for the high T in master
        ψ_ref, E_ref[t₀+1:T+1], sim_ref, mcs_ref = nMCS_energyDynamic(ψ_ref, sim_ref, T-t₀, ADJUST_INTERVAL)
        
        #Gather results from workers
        for i=1:NWS
            ψ_w[i], E_w[i,t₀+1:T+1], sim_w[i], mcs_w[i] = fetch(ψ_future_list[i])
        end
        adjustment_mcs += max(maximum(mcs_w),mcs_ref)
        tₛ=T
        
        #Check if for all processes dE < 0, and get largest step that it happens for
        for i=t₀+1:T+1
            for w=1:3
                if ( E_ref[i] - E_w[w,i] <= 0.0 && E_check_workers[w] == 0 )
                    E_check_workers[w] = i
                    println("Worker $(w) initially thermalised after $(i) steps")
                    flush(STDOUT)
                end
            end
        end

        if minimum(E_check_workers) == 0 #This is always true if we found no dE < 0 in all workers
            t₀ = T
            T = min(Int(ceil(T*ex)), CUTOFF_MAX)
            println("Increasing simulation time with t0 = $(t₀) and T = $(T)")
            if T == CUTOFF_MAX
                println("Failed to find a point where ΔE <= 0")
                return (-1, E_ref[1:T], E_w[:,1:T], ψ_ref, ψ_w, sim_ref, sim_w)
            end
        else
            t₀ = maximum(E_check_workers)
            println("All workers initially thermalised after $(t₀) steps")
            thermalised_init = true
        end
    end

    tₛ = T
    T = T + AVG_MCS_INTERVAL
    #Set this to A_I so that the simulation constants are only updated by energyDynamic at the END of each 
    #iteration in the while loop.
    #Increase simulation time by ADJUST_INT_THERM and try to find an interval with small average dE
    while T<CUTOFF_MAX
        #Start by updating the simulation constants
        for i=1:NWS
            ψ_future_list[i] = @spawn(adjustSimConstantsPar(sim_w[i], ψ_w[i]))
        end
        ψ_ref, sim_ref, mcs_ref = adjustSimConstantsPar(sim_ref, ψ_ref)
        for i=1:NWS
            ψ_w[i], sim_w[i], mcs_w[i] = fetch(ψ_future_list[i]) 
        end
        adjustment_mcs += max(maximum(mcs_w),mcs_ref)
        
        #Do MCS to find average
        for i=1:NWS
            ψ_future_list[i] = @spawn nMCS_energy(ψ_w[i], sim_w[i], AVG_MCS_INTERVAL)
        end
        ψ_ref, E_ref[tₛ+1:T+1] = nMCS_energy(ψ_ref, sim_ref, AVG_MCS_INTERVAL)
        for i=1:NWS
            ψ_w[i], E_w[i,tₛ+1:T+1] = fetch(ψ_future_list[i])
        end

        if check_thermalisation(E_ref/N, E_w/N, NWS, tₛ, T,STD_NUMBER)
            println("Thermalisation succesfull")
            return(T+adjustment_mcs, E_ref[1:T], E_w[:,1:T], ψ_ref, ψ_w, sim_ref, sim_w)
        end
        tₛ=T
        T = T + AVG_MCS_INTERVAL
        println("Failed to find proper therm, increasing simulation time to T = $(T)")
    end
    return -1, E_ref, E_w, ψ_ref, ψ_w, sim_ref, sim_w
end
