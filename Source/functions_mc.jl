####################################################################################################
#                            Monte-Carlo functions
#
####################################################################################################

# --------------------------------------------------------------------------------------------------
# Given a lattice site ϕ, propose a new lattice site with values in intervals around the existing ones.
function proposeLocalUpdate(ϕ::LatticeSite, sim::Controls)
    UMAX::Int64 = 4
    u⁺ = mod(ϕ.u⁺ + rand(Uniform(-sim.umax,sim.umax)), UMAX) # This does not allow u⁺ = UMAX, is this a problem?
	u⁻ = 0.0#mod(ϕ.u⁻ + rand(Uniform(-sim.umax,sim.umax)), UMAX)
    # Construct new configuration at lattice site.
    #return LatticeSite([ϕ.A[1]+rand(Uniform(-sim.Amax,sim.Amax)), ϕ.A[2]+rand(Uniform(-sim.Amax,sim.Amax)),
    #                    ϕ.A[3]+rand(Uniform(-sim.Amax,sim.Amax))],
    #    mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
    #    u⁺, u⁻)
#    return LatticeSite([ϕ.A[1]+rand(Uniform(-sim.Amax,sim.Amax)), ϕ.A[2]+rand(Uniform(-sim.Amax,sim.Amax)),
#                        ϕ.A[3]+rand(Uniform(-sim.Amax,sim.Amax))],
#        mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
#        u⁺, u⁻)
    #return LatticeSite([0, 0, 0],
    #    mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
    #    u⁺, u⁻)
    return LatticeSite([0, 0, 0],
        mod(ϕ.θ⁺ + rand(sim.θ_rng), 2π), mod(ϕ.θ⁻ + rand(sim.θ_rng), 2π), 
        u⁺, u⁻)
end

# --------------------------------------------------------------------------------------------------
# Performes a Metropolis Hasting update on a lattice site at position pos in state ψ given an inverse temperature
# β and where ϕᵣ... gives nearest and next nearest neighbor sites. Note that pos gives [y,x] of the position of
# the lattice site in normal array notation such that [1,1] is the upper left corner.
function metropolisHastingUpdate!(ψ::State, pos::Array{Int64,1}, sim::Controls)
	# Save the lattice site at the targeted position in a temporary variable ϕ and use the lattice site
	# as a basis for proposing a new lattice site ϕ′. Then find the energy difference between having
	# ϕ′ or ϕ at position pos.
    ϕ = ψ.lattice[pos...]
    ϕ′ = proposeLocalUpdate(ϕ, sim)
    δE = ΔE(ϕ′, ϕ, ψ.nb[pos...], ψ.nnb[pos...], ψ.nnnb[pos...], pos[2], ψ.consts)
    
    # Create random number ran ∈ (0,1].
    ran = 1-rand()
#    if ran==0.0
#        ran=1.0
#    end
    
    # Update state with probability min(1, e^{-β⋅δE})
    # and return the energy of final state regardless of whether it gets updated or not.
    if log(ran) <= -ψ.consts.β*δE
		set!(ϕ, ϕ′)
        return δE
    else
        return 0.0
    end
end

# --------------------------------------------------------------------------------------------------
# Takes a state ψ with an L×L lattice and tries to update each site on the lattice by running the
# metropolisHastingUpdate! function on it. Each part of the boundary is updated separately so that periodic
# boundary conditions are taken care of for values stored in each lattice site.
function mcSweep!(ψ::State, sim::Controls = Controls(π/3, 0.4, 3.0))
   
    # Find size of the lattice L
    L::Int64 = ψ.consts.L
    L₃::Int64 = ψ.consts.L₃
    
    for z_pos = 1:L₃, h_pos = 1:L, v_pos = 1:L
        metropolisHastingUpdate!(ψ, [v_pos,h_pos,z_pos], sim)
    end
end

# --------------------------------------------------------------------------------------------------
# Same as above but returns the fraction of accepted proposals over number of proposals and calculates
# nearest neighbors in a dynamic but more costly way.
function mcSweepFrac!(ψ::State, sim::Controls = Controls(π/3, 0.4, 1.0))
    count = 0
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    for z_pos=1:L₃, h_pos=1:L, v_pos=1:L
        if metropolisHastingUpdate!(ψ,[v_pos,h_pos,z_pos], sim) != 0.0
            count += 1
        end
    end

    return count/(L^2*L₃)
end


# -----------------------------------------------------------------------------------------------------------
# Same as mcSweep but adds up the energy-differences and returns it.
function mcSweepEn!(ψ::State, sim::Controls = Controls(π/3, 0.4, 3.0))
    
    δE = 0.0
    # Find size of the lattice L
    L::Int64 = ψ.consts.L
    
    # Update the bulk
    for z_pos=1:ψ.consts.L₃, h_pos=1:L, v_pos=1:L
        δE += metropolisHastingUpdate!(ψ, [v_pos,h_pos,z_pos], sim)
    end
    return δE
end


# -------------------------------------------------------------------------------------------------
# Takes a state ψ and plots the fraction of accepted proposals on number of proposals for M
# Monte Carlo sweeps of the lattice. Then return the average and standard deviation of this
# fraction as well as the time-series itself.
function mcProposalFraction(ψ::State, sim::Controls=Controls(π/2, 0.4, 3.0), M::Int64=500)
	ψ_copy = copy(ψ)
    L = ψ.consts.L
    fracs = zeros(M)
    
    # Go through the entire lattice M times and gain the statistic of whether it gets updated or not
    for i = 1:M
        fracs[i] = mcSweepFrac!(ψ_copy, sim)
    end
    res = mean(fracs)
    stdev = std(fracs)
    return (res, stdev, fracs)
end

# Same as mcProposalFraction but does update the state.
# The number of times the state has been updated can be surmised from the input parameter M.
# Does not return the series of acceptance-rates unlike mcProposalFraction
function mcProposalFraction!(ψ::State, sim::Controls=Controls(π/2, 0.4, 3.0), M::Int64=500)
    L = ψ.consts.L
    fracs = zeros(M)
    
    # Go through the entire lattice M times and gain the statistic of whether it gets updated or not
    for i = 1:M
        fracs[i] = mcSweepFrac!(ψ, sim)
    end
    av = mean(fracs)
    stdev = std(fracs)
    return (av, stdev)
end

# --------------------------------------------------------------------------------------------------
function nMCS(ψ::State, sim::Controls, n::Int64)
    for i=1:n
        mcSweep!(ψ, sim)
    end
    return ψ
end

# --------------------------------------------------------------------------------------------------
# Preform nMCS for all states in a list, as much as possible in parallell.
function nMCS!(ψ_list::Array{State, 1}, sim_list::Array{Controls, 1}, N::Int64)
    nw = nprocs()-1
    n_state = length(ψ_list)
    
    i = 0 # Index in ψ_list of states already updated.
    while i < n_state
        
        worker_jobs = min(nw, n_state-1-i) # Number of needed jobs given to workers
        # Start the max number of workers if that wouldn't be too much.
        work_futures = [Future() for w = 1:worker_jobs]
        
        for w = 1:worker_jobs
            work_futures[w] = @spawn nMCS(ψ_list[i+w], sim_list[i+w], N)
        end
        nMCS(ψ_list[i+worker_jobs+1], sim_list[i+worker_jobs+1], N)
        
        for w = 1:worker_jobs
            ψ_list[i+w] = fetch(work_futures[w])
        end
        
        i += worker_jobs+1
    end
    
    # After this, all states should have been updated
    return ψ_list
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

# --------------------------------------------------------------------------------------------------
# Preform nMCSEnergy on a list of states (as much as possible in paralllel)
function nMCSEnergy!(ψ_list::Array{State, 1}, sim_list::Array{Controls, 1}, n::Int64, E₀_list::Array{Float64, 1})
    nw = nprocs()-1
    n_state = length(ψ_list)
    E_matrix = Array{Float64, 2}(n_state, n)
    
    i = 0 # Index in ψ_list of states already updated.
    while i < n_state
        
        worker_jobs = min(nw, n_state-1-i) # Number of needed jobs given to workers
        # Start the max number of workers if that wouldn't be too much.
        work_futures = [Future() for w = 1:worker_jobs]
        
        for w = 1:worker_jobs
            work_futures[w] = @spawn nMCSEnergy(ψ_list[i+w], sim_list[i+w], n, E₀_list[i+w])
        end
        index = i+worker_jobs+1
        ψ_list[index], E_matrix[index, :] =  nMCSEnergy(ψ_list[index], sim_list[index], n, E₀_list[index])
        
        for w = 1:worker_jobs
            ψ_list[i+w], E_matrix[i+w, :] = fetch(work_futures[w])
        end
        
        i += worker_jobs+1
    end
    
    # After this, all states should have been updated
    return ψ_list, E_matrix
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



# -----------------------------------------------------------------------------------------------------------
# Adjusts the simulation controls such that ψ has an AR larger than LOWER. Tries to find the simulation controls
# that are simultaneously the ones with the highest values.
function adjustSimConstants!(sim::Controls, ψ::State, M::Int64 = 40)
    #println("Adjusting simulation constants $(sim.θmax), $(sim.umax), $(sim.Amax)")
    CUTOFF_MAX = 42       # How many times the while loop should run. 
    LOWER = 0.3           # Minimum acceptance rate (AR)
    NEEDED_PROPOSALS = 5  # Number of sim with AR >= LOWER
    DIVI = 1.5            # The number the current value is divided by to get lower limit of search interval.
    TRIED_VALUES = 4      # Number of values to try in new interval
    
    proposedConstants = [copy(sim) for i=1:NEEDED_PROPOSALS]
    proposedAR = zeros(NEEDED_PROPOSALS)
    proposals = 0
    n = 0
    adjustment_mcs = 0
    
    # First we get an estimate of the accept probability
    (av, std) = mcProposalFraction!(ψ, sim, M)
    adjustment_mcs = M
    if av >= LOWER
        return (av, adjustment_mcs, ψ, sim)
    end
    
    s₀ = copy(sim)
    while (proposals < NEEDED_PROPOSALS)
        # The starting state of this loop will be that we have no proposals for sims that has acceptance rate higher
        # than LOWER (including the initial s₀)
        
        # First we look at Amax
        interval_end = s₀.Amax/DIVI
        tries_A = [Controls(s₀.θmax, s₀.umax, 
                s₀.Amax-x*(s₀.Amax-interval_end)/TRIED_VALUES) for x = 1:TRIED_VALUES]
        f_A = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_A[i], std) = mcProposalFraction!(ψ, tries_A[i], M)
            adjustment_mcs += M
            # If we find a value with probability >= LOWER, then this is the largest such value we have found
            # and should be included in proposed constants
            if f_A[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_A[i]
                proposedAR[proposals] = f_A[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break

        # Then we try to vary umax
        interval_end = s₀.umax/DIVI
        tries_u = [Controls(s₀.θmax, s₀.umax - x*(s₀.umax-interval_end)/TRIED_VALUES,
                s₀.Amax) for x = 1:TRIED_VALUES]
        f_u = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_u[i], std) = mcProposalFraction!(ψ, tries_u[i], M)
            adjustment_mcs += M
            if f_u[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_u[i]
                proposedAR[proposals] = f_u[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break

        # Then we try to vary θmax
        interval_end = s₀.θmax/DIVI
        tries_θ = [Controls(s₀.θmax - x*(s₀.θmax - interval_end)/TRIED_VALUES, s₀.umax,
                s₀.Amax) for x = 1:TRIED_VALUES]
        f_θ = zeros(TRIED_VALUES)
        for i = 1:TRIED_VALUES
            (f_θ[i], std) = mcProposalFraction!(ψ, tries_θ[i], M)
            adjustment_mcs += M
            if f_θ[i] >= LOWER
                proposals += 1
                proposedConstants[proposals] = tries_θ[i]
                proposedAR[proposals] = f_θ[i]
                break
            end
        end
        
        # If we have the needed number of proposals we exit the loop.
        (proposals >= NEEDED_PROPOSALS) && break
        
        # At this point in the loop we have tried to get some proposals but failed to get enough of them.
        # We need to start the loop again with a new starting state that is such that the accept ratio is
        # as high is possible so that we can get more proposals.
        accept_ratios = vcat(f_A, f_u, f_θ)
        tries = vcat(tries_A, tries_u, tries_θ)
        i_max = indmax(accept_ratios)   # Finding index of sim that gave highest accept ratio.
        s₀ = tries[i_max]               # Setting this sim to the initial one.
        
        
        n += 1
        if n >= CUTOFF_MAX
            println("WARNING: Could not find simulation constant such that update probability 
                was higher than $(LOWER)")
            sim = s₀
            return (accept_ratios[i_max], adjustment_mcs, ψ, sim)
        end
    end
    # The end situation of the loop is that we have a number of proposals >= NEEDED_PROPOSALS.
    
    # Finding the distance from zero of the different proposals.
    norms = zeros(proposals)
    for i = 1:proposals
        norms[i] = proposedConstants[i].θmax^2 + proposedConstants[i].umax^2 + proposedConstants[i].Amax^2
    end
    i_max = indmax(norms)
    setValues!(sim, proposedConstants[i_max]) # Finally we update the simulation constants to the sim that has 
    # highest norm, and an accept ratio above LOWER.
    return (proposedAR[i_max], adjustment_mcs, ψ, sim)
    # Return the acceptance ratio of the new sim and the number of Monte-Carlo Sweeps done during this adjustment.
end

# --------------------------------------------------------------------------------------------------
# Adjust controls of all states in a list (as much as possible in parallel)
function adjustSimConstants!(ψ_list::Array{State, 1}, sim_list::Array{Controls, 1})
    nw = nprocs()-1
    n_state = length(ψ_list)
    mcs_list = zeros(Int64, n_state)
    ar_list = Array{Float64, 1}(n_state)
    
    i = 0 # Index in ψ_list of states already updated.
    while i < n_state
        
        worker_jobs = min(nw, n_state-1-i) # Number of needed jobs given to workers
        # Start the max number of workers if that wouldn't be too much.
        work_futures = [Future() for w = 1:worker_jobs]
        
        for w = 1:worker_jobs
            work_futures[w] = @spawn adjustSimConstants!(sim_list[i+w], ψ_list[i+w])
        end
        index = i+worker_jobs+1
        ar_list[index], mcs_list[index], ψ_list[index], sim_list[index] = adjustSimConstants!(sim_list[index],
            ψ_list[index])
        
        # Fetch futures
        for w = 1:worker_jobs
            ar_list[i+w], mcs_list[i+w], ψ_list[i+w], sim_list[i+w] = fetch(work_futures[w])
        end
        
        i += worker_jobs+1
    end
    
    # After this, all states should have been updated
    return mcs_list, ar_list
end
