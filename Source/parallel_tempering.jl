# Implementation of parallel tempering with possibility of extension to FOPT through histograms.
# Assumes that mcSweepEn! and DList is available.
# Use: Given a set of states in ψ_list with different temperatures, but otherwise same system constants, construct a PTRun
# with pt = PTRun(ψ_list, sim_list, N_mc). N_mc is normally 1. A parallel tempering step is then preformed with
# PTStep!(pt). The up-mover histograms can be obtained with getHistograms(pt), and the acceptance rates for global updates
# between neighboring temperatures can be gotten with getARList(pt).

@everywhere mutable struct Replica
    ψ::State
    state::Int64
    En::Float64
    sim::Controls
end

function Replica(ψ::State, sim::Controls)
    if ψ.PBC
        out = Replica(ψ, 3, E(ψ), sim)
    else
        out = Replica(ψ, 3, E_APBC(ψ), sim)
    end
    return out
end

mutable struct PTRun
    dr_list::DList              # A list of replicas distributed over available processes.
    accepts::Array{Int64,1}     # Number of accepted PT-swaps between i and i+1
    histograms::Array{Array{Int64, 1},1}
    # This is a list of histograms; 3 histograms pr. temperature. The first is the number of
    # up-moving visits, the 2nd is the number of down-moving visits and the third is the number of visits
    # of uninitialized states.
    β_list::Array{Float64, 1}   # List of inverse temperatures
    rep_map::Array{Int64, 1}       # A map from the index of the inverse temp. in β_list to the corresponding index in
    # dr_list that has the replica with the corresponding temperature.
    E_list::Array{Float64, 1}   # List of energies of the replicas. The index in this list should mirror the indices
                                # in dr_list.
    N_mc::Int64                 # Number of MC sweeps to do between each PT step.
    N_pt::Int64                 # Number of PT steps done
    N_temp::Int64               # Number of temperatures.
end

function PTRun(ψ_list::Array{State,1}, sim_list::Array{Controls,1}, N_mc::Int64; verbose=false)
    N = length(ψ_list)
    N == length(sim_list) || throw(error("ERROR: List of states and control constants are not same length"))
    
    # Sort ψ_list according to increasing temperature
    β_list = [ψ.consts.β for ψ in ψ_list]
    perm = sortperm(β_list; rev=true)
    if verbose && perm != [i for i = 1:N]
        println("WARNING: Inserted state list did not have states in increasing temperature")
    end
    β_list = β_list[perm]
    ψ_list = ψ_list[perm]
    sim_list = sim_list[perm]
    
    # Initialize remaining lists.
    accepts = [0 for i = 1:N-1]
    rep_list = [Replica(ψ_list[i], sim_list[i]) for i = 1:N]
    rep_list[1].state = 1    # The lowest temperature replica is moving up
    rep_list[N].state = 2    # The highest temperature replica is moving down.
    rep_map = [i for i = 1:N]
    E_list = [R.En for R in rep_list]
    histograms = [[0,0,0] for i = 1:N]
    dr_list = DList(rep_list)
    return PTRun(dr_list, accepts, histograms, β_list, rep_map, E_list, N_mc, 0, N)
end

function E(pt::PTRun)
    return pt.E_list[pt.rep_map]
end

@everywhere function nMCS!(R::Replica, n::Int64)
    En = R.En
    for i=1:n
        En += mcSweepEn!(R.ψ,R.sim)
    end
    R.En = En
    return R
end

@everywhere function distributeTemperatures!(chan::RemoteChannel{Channel{Array{Replica, 1}}}, β_list::Array{Float64, 1})
    local_list = take!(chan)
    for (i, β) = enumerate(β_list)
        s = local_list[i].ψ.consts
        new_syst = SystConstants(s.L, s.L₃, s.g⁻², s.ν, s.νₘ, s.νₐₙ, s.κ₅, s.f, β)
        local_list[i].ψ.consts = new_syst
    end
    put!(chan, local_list)
    nothing
end
function distributeTemperatures!(dr_list::DList, β_list::Array{Float64, 1})
    futures = Array{Future, 1}(undef, dr_list.chuncks)
    for (i, ck) = enumerate(dr_list.chunck_list)
        futures[i] = @spawnat ck.p distributeTemperatures!(ck.chan, β_list[ck.range])
    end
    for i = 1:dr_list.chuncks
        wait(futures[i])
    end
    nothing
end

function distributeTemperatures!(pt::PTRun, β_list::Array{Float64, 1})
    pt.β_list .= β_list
    distributeTemperatures!(pt.dr_list, pt.β_list)
    nothing
end

function incrementHistograms!(pt::PTRun)
    states = dmap(R -> R.state, pt.dr_list)
    for i = 1:pt.N_temp
        histogram = pt.histograms[i]
        replica_state = states[pt.rep_map[i]]
        histogram[replica_state] += 1
    end
end

function dmap(f::Function, pt::PTRun)
    return dmap(f, pt.dr_list)[pt.rep_map]
end

function ptSwap!(pt::PTRun, i::Int64, j::Int64)
    # i and j refer to the indices in the temperature array. We swap the contents of rep_map at these indices.
    # This means that effectively we swap the states associated with the temperatures at i and j.
    temp = pt.rep_map[j]
    pt.rep_map[j] = pt.rep_map[i]
    pt.rep_map[i] = temp
    
    return nothing
end

function attemptSwap(pt::PTRun, i::Int64, j::Int64)
    
    Δβ = pt.β_list[j] - pt.β_list[i]#ψⱼ.consts.β-ψᵢ.consts.β
    ΔE = pt.E_list[pt.rep_map[j]] - pt.E_list[pt.rep_map[i]]#rep_list[j].En - rep_list[i].En
    
    r = 1-rand() # make a number ∈ (0, 1]
    # Same method as metropolis-hasting selection: to select with probability min{1, e^{Δβ⋅ΔE}}
    # We could choose a random number r ∈ [0, 1] and check if it's lower than r <= e^{Δβ⋅ΔE} provided
    # the exponent is negative. Taking the logarithm on both sides we see that in this case then also
    # log(r) <= Δβ⋅ΔE. This has the advantage that since r ∈ (0,1], log(r) will always be negative
    # and thus we always select if Δβ⋅ΔE is positive, as we should.
    if log(r) <= Δβ*ΔE
        ptSwap!(pt, i, j)
        return true
    end
    return false
end

@everywhere function setState(chan::RemoteChannel{Channel{Array{Replica, 1}}}, i::Int64, state::Int64)
    rep_list = take!(chan)
    rep_list[i].state = state
    put!(chan, rep_list)
    nothing
end
function setStates(dr_list::DList, up_index::Int64, down_index::Int64)
    up_tup = dr_list.chunck_map[up_index]
    up_ck = dr_list.chunck_list[up_tup[1]]
    up_fut = @spawnat up_ck.p setState(up_ck.chan, up_tup[2], 1)
    
    down_tup = dr_list.chunck_map[down_index]
    down_ck = dr_list.chunck_list[down_tup[1]]
    down_fut = @spawnat down_ck.p setState(down_ck.chan, down_tup[2], 2)
    
    wait(up_fut)
    wait(down_fut)
    nothing
end

function PTStep!(pt::PTRun)
    # First we preform N_mc MC sweeps for all of the replicas
    N_mc = pt.N_mc
    dmutate(R -> nMCS!(R, N_mc), pt.dr_list)
    # Get energies
    pt.E_list = dmap(R -> R.En, pt.dr_list)
    
    # Before a squence of swap moves we increment the histograms
    incrementHistograms!(pt)
    
    # Then we go through the temperatures and make a swap with probability
    # p = min{1, e^(-Δβ*ΔE)} between neighboring temperatures.
    for i = 1:pt.N_temp-1
        if attemptSwap(pt, i, i+1)
            # Updates acceptance rates if successful
            pt.accepts[i] += 1
        end
    end
    
    # Distribute the new temperatures to worker processes
    distributeTemperatures!(pt.dr_list, pt.β_list[pt.rep_map[pt.rep_map]])
    
    # Finally increment number of PT steps done and set the end states
    # Since the lowest temperature is at index 1, we set the replica at rep_map[1] to be a up-mover.
    # The temperature at N_temp is highest, thus the replica at rep_map[N_temp] is set to a down-mover.
    pt.N_pt += 1
    setStates(pt.dr_list, pt.rep_map[1], pt.rep_map[pt.N_temp])
    
    nothing
end

function getARList(pt::PTRun)
    return [an/(pt.N_pt) for an in pt.accepts]
end

function getHistograms(pt::PTRun)
    return [hist[1]/(hist[1]+hist[2]) for hist in pt.histograms]
end

