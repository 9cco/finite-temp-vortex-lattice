############################################################################################################################
#                               Metropolis-Hastings update of LatticeSites // SubCuboidModule
#__________________________________________________________________________________________________________________________#
############################################################################################################################


function proposeLocalUpdate(ϕ::LatticeSite, sim::Controls)
    UMAX::Int64 = 4
    u⁺ = 1/√(2)#√(0.2)#mod(ϕ.u⁺ + rand(Uniform(-sim.umax,sim.umax)), UMAX) # This does not allow u⁺ = UMAX, is this a problem?
    u⁻ = 1/√(2)#√(2)#mod(ϕ.u⁻ + rand(Uniform(-sim.umax,sim.umax)), UMAX)
    # Construct new configuration at lattice site.
    #return LatticeSite([ϕ.A₁+rand(Uniform(-sim.Amax,sim.Amax)), ϕ.A₂+rand(Uniform(-sim.Amax,sim.Amax)),
    #                    ϕ.A₃+rand(Uniform(-sim.Amax,sim.Amax))],
    #    mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
    #    u⁺, u⁻)
    return LatticeSite(ϕ.A₁+rand(sim.A_rng), ϕ.A₂+rand(sim.A_rng),
                        ϕ.A₃+rand(sim.A_rng),
        (2*rand()-1)*π, (2*rand()-1)*π, # Get θ ∈ [-π, π)
    #   mod(ϕ.θ⁺ + rand(sim.θ_rng), 2π)-π, mod(ϕ.θ⁻ + rand(sim.θ_rng), 2π)-π, #0.0,#
        u⁺, u⁻, ϕ.x)
    #return LatticeSite([0, 0, 0],
    #    mod(ϕ.θ⁺ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), mod(ϕ.θ⁻ + rand(Uniform(-sim.θmax,sim.θmax)), 2π), 
    #    u⁺, u⁻)
    #return LatticeSite([0.0, 0.0, 0.0],
    #    mod(ϕ.θ⁺ + rand(sim.θ_rng), 2π), mod(ϕ.θ⁻ + rand(sim.θ_rng), 2π), 
    #    u⁺, u⁻)
end

# A dummy function that should be transplanted with proper Metropolis-Hastings update in final version
#@everywhere function updateLatticeSite!(ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors, β::R) where R<:Real
#    ϕ.u⁺ = β
#    true
#end

# The real function for updating LatticeSites
function updateLatticeSite!(ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors,
        syst::SystConstants, sim::Controls, β::R) where R<:Real
    
    ϕ′ = proposeLocalUpdate(ϕ, sim)
    δE = ΔE(ϕ′, ϕ, nb, nnb, syst)
    
    # Update following a Metrolopis-Hastings selection.
    ran = 1 - rand()
    
    if log(ran) <= -β*δE
        set!(ϕ, ϕ′)
        return true
    else
        return false
    end
    #set!(ϕ,ϕ′)
    #return (true, 1.0)
end
function updateLatticeSiteRetEn!(ϕ::LatticeSite, nb::NearestNeighbors, nnb::NextNeighbors,
        syst::SystConstants, sim::Controls, β::R) where R<:Real
    
    ϕ′ = proposeLocalUpdate(ϕ, sim)
    δE = ΔE(ϕ′, ϕ, nb, nnb, syst)
    
    # Update following a Metrolopis-Hastings selection.
    ran = 1 - rand()
    
    if log(ran) <= -β*δE
        set!(ϕ, ϕ′)
        return true, δE
    else
        return false, 0.0
    end
    #set!(ϕ,ϕ′)
    #return (true, 1.0)
end


