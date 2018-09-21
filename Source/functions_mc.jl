####################################################################################################
#                            Monte-Carlo functions
#
####################################################################################################

# --------------------------------------------------------------------------------------------------
# Given a lattice site ϕ, propose a new lattice site with values in intervals around the existing ones.
function proposeLocalUpdate(ϕ::LatticeSite, sim::Controls)
    
    u⁺ = mod(ϕ.u⁺ + rand(Uniform(-sim.umax,sim.umax)),1) # This does not allow u⁺ = 1, is this a problem?
    # Construct new configuration at lattice site.
    return LatticeSite([ϕ.A[1]+rand(Uniform(-sim.Amax,sim.Amax)), ϕ.A[2]+rand(Uniform(-sim.Amax,sim.Amax))],
        mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
        u⁺, √(1-u⁺^2))
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
    ran = rand()
    if ran==0
        ran=1
    end
    
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
    
    for h_pos = 1:L, v_pos = 1:L
        metropolisHastingUpdate!(ψ, [v_pos,h_pos], sim)
    end
end

# --------------------------------------------------------------------------------------------------
# Same as above but returns the fraction of accepted proposals over number of proposals and calculates
# nearest neighbors in a dynamic but more costly way.
function mcSweepFrac!(ψ::State, sim::Controls = Controls(π/3, 0.4, 1.0))
    count = 0
    L = ψ.consts.L
    for h_pos = 1:L
        for v_pos = 1:L
            if metropolisHastingUpdate!(ψ,[v_pos,h_pos], sim) != 0.0
                count += 1
            end
        end
    end
    return count/L^2
end


# -----------------------------------------------------------------------------------------------------------
# Same as mcSweep but adds up the energy-differences and returns it.
function mcSweepEn!(ψ::State, sim::Controls = Controls(π/3, 0.4, 3.0))
    
    δE = 0.0
    # Find size of the lattice L
    L::Int64 = ψ.consts.L
    
    # Update the bulk
    for h_pos=1:L, v_pos=1:L
        δE += metropolisHastingUpdate!(ψ, [v_pos,h_pos], sim)
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
        return (av, adjustment_mcs)
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
            return (accept_ratios[i_max], adjustment_mcs)
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
    return (proposedAR[i_max], adjustment_mcs)
    # Return the acceptance ratio of the new sim and the number of Monte-Carlo Sweeps done during this adjustment.
end


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



