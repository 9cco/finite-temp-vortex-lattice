using Distributed
using Distributions
using Test
using BenchmarkTools
#using Plots
#pyplot()

############################################################################################################################
#                               Types
#__________________________________________________________________________________________________________________________#
############################################################################################################################

const two_pi = 2π
struct CubConstants
    l₁::Int64; l₂::Int64; l₃::Int64
end
struct SystConstants
    L₁::Int64     # System size in x-direction
    L₂::Int64     # System size in y-direction
    L₃::Int64     # System size in z-direction
    g⁻²::Float64  # 1/g² for gauge coupling g
    ν::Float64    # Anisotropy constant
    κ₅::Float64   # z-layer coupling.
    f::Float64    # Magnetic filling fraction
end
mutable struct LatticeSite
    A₁::Float64  # Integrated fluctuating vector potential x-direction
    A₂::Float64 # y-direction
    A₃::Float64 # z-direction
    θ⁺::Float64 # Phase of the + component
    θ⁻::Float64 # Phase of the - component
    u⁺::Float64 # Amplitude of + component
    u⁻::Float64 # Amplitude of - component (should always be √u⁺)
    x::Int64   # Position in x-direction
end
mutable struct Controls
    θmax::Float64
    umax::Float64
    Amax::Float64
    θ_rng::Distributions.Uniform{Float64}
    u_rng::Distributions.Uniform{Float64}
    A_rng::Distributions.Uniform{Float64}
end
struct NearestNeighbors
	ϕᵣ₊₁::LatticeSite
	ϕᵣ₋₁::LatticeSite
	ϕᵣ₊₂::LatticeSite
	ϕᵣ₋₂::LatticeSite
    ϕᵣ₊₃::LatticeSite
    ϕᵣ₋₃::LatticeSite
end
struct NextNeighbors
	ϕᵣ₊₁₋₂::LatticeSite
	ϕᵣ₋₁₊₂::LatticeSite
    ϕᵣ₊₁₋₃::LatticeSite
    ϕᵣ₋₁₊₃::LatticeSite
    ϕᵣ₊₂₋₃::LatticeSite
    ϕᵣ₋₂₊₃::LatticeSite
end
struct RemoteNeighbors{T}
    rnᵣ₊₁::RemoteChannel{Channel{T}}
    rnᵣ₋₁::RemoteChannel{Channel{T}}
    rnᵣ₊₂::RemoteChannel{Channel{T}}
    rnᵣ₋₂::RemoteChannel{Channel{T}}
    rnᵣ₊₃::RemoteChannel{Channel{T}}
    rnᵣ₋₃::RemoteChannel{Channel{T}}
    rnᵣ₊₁₋₂::RemoteChannel{Channel{T}}
    rnᵣ₊₁₋₃::RemoteChannel{Channel{T}}
    rnᵣ₋₁₊₂::RemoteChannel{Channel{T}}
    rnᵣ₋₁₊₃::RemoteChannel{Channel{T}}
    rnᵣ₊₂₋₃::RemoteChannel{Channel{T}}
    rnᵣ₋₂₊₃::RemoteChannel{Channel{T}}
end
# A sub-cuboid lives on an (external) process. It contains a shell of nearest neighbor lattice sites
# in addition to the lattice sites that it contains.
mutable struct SubCuboid
    lattice::Array{LatticeSite, 3}    # The lattice sites of the cuboid
    nb::Array{NearestNeighbors, 3}    # References to nearest neighbors of the lattice sites in lattice
    nnb::Array{NextNeighbors, 3}      # References to upper-left, lower-right diagonal neighbors
    sc_nb::RemoteNeighbors{SubCuboid}      # Each sub-cuboid should have a RemoteReference to all the neighboring
                                      # SubCuboids so that updated values can be sent directly.
    shell::Vector{Array{LatticeSite, 2}}      # Shell around the sub cuboid consisting of nearest neighbors of the lattice
                                      # sites on the perifery of the cuboid. These have equal values to lattice sites
                                      # at neighboring sub-cuboids. Contains 6 planes of neighbor points.
    shell_edge14::Vector{LatticeSite} # Two vectors of lattice sites for the 41 and 23 edge.
    shell_edge16::Vector{LatticeSite}
    shell_edge23::Vector{LatticeSite}
    shell_edge25::Vector{LatticeSite}
    shell_edge36::Vector{LatticeSite}
    shell_edge45::Vector{LatticeSite}
    consts::CubConstants
    syst::SystConstants
    sim::Controls
    β::Float64
end


############################################################################################################################
#                            Functions for ::Controls
#__________________________________________________________________________________________________________________________#
############################################################################################################################

function Controls(θ_max::Float64, u_max::Float64, A_max::Float64)
    return Controls(θ_max, u_max, A_max, Uniform(-θ_max,θ_max), Uniform(-u_max,u_max), Uniform(-A_max, A_max))
end

function Controls()
    return Controls(π/3, 0.4, 3.0)
end

import Base.copy
function copy(sim::Controls)
    return Controls(sim.θmax, sim.umax, sim.Amax)
end

# -------------------------------------------------------------------------------------------------
function setValues!(sim_target::Controls, sim_source::Controls)
    sim_target.θmax = sim_source.θmax
    sim_target.umax = sim_source.umax
    sim_target.Amax = sim_source.Amax
    sim_target.θ_rng = Uniform(-sim_source.θmax,sim_source.θmax)
    sim_target.u_rng = Uniform(-sim_source.umax,sim_source.umax)
    sim_target.A_rng = Uniform(-sim_source.Amax,sim_source.Amax)
    return 
end

############################################################################################################################
#                               Neighbor functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

include("neighbors.jl")

############################################################################################################################
#                               Utility functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# -------------------------------------------------------------------------------------------------
# LatticeSite outer constructor
# Initializes a lattice site that has radom values at a horizontal position.
function LatticeSite(h_pos::Int64)
    A_max = 3.0
    u⁺ = rand()
    LatticeSite(rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max,A_max)),
                two_pi*rand(), two_pi*rand(), u⁺, √(1-u⁺^2), h_pos)
end

import Base.copy
function copy(re_nb::RemoteNeighbors{T}) where T
    rnᵣ₊₁ = re_nb.rnᵣ₊₁
    rnᵣ₋₁ = re_nb.rnᵣ₋₁
    rnᵣ₊₂ = re_nb.rnᵣ₊₂
    rnᵣ₋₂ = re_nb.rnᵣ₋₂
    rnᵣ₊₃ = re_nb.rnᵣ₊₃
    rnᵣ₋₃ = re_nb.rnᵣ₋₃
    rnᵣ₊₁₋₂ = re_nb.rnᵣ₊₁₋₂
    rnᵣ₊₁₋₃ = re_nb.rnᵣ₊₁₋₃
    rnᵣ₋₁₊₂ = re_nb.rnᵣ₋₁₊₂
    rnᵣ₋₁₊₃ = re_nb.rnᵣ₋₁₊₃
    rnᵣ₊₂₋₃ = re_nb.rnᵣ₊₂₋₃
    rnᵣ₋₂₊₃ = re_nb.rnᵣ₋₂₊₃
    RemoteNeighbors{T}(rnᵣ₊₁, rnᵣ₋₁, rnᵣ₊₂, rnᵣ₋₂, rnᵣ₊₃, rnᵣ₋₃, rnᵣ₊₁₋₂, rnᵣ₊₁₋₃, rnᵣ₋₁₊₂, rnᵣ₋₁₊₃, rnᵣ₊₂₋₃, rnᵣ₋₂₊₃)
end

# -------------------------------------------------------------------------------------------------
# Mutates target to have the same values as src.
function set!(target::LatticeSite, src::LatticeSite)
    target.A₁ = src.A₁
    target.A₂ = src.A₂
    target.A₃ = src.A₃
    target.θ⁺ = src.θ⁺
    target.θ⁻ = src.θ⁻
    target.u⁺ = src.u⁺
    target.u⁻ = src.u⁻
    target.x = src.x
    nothing
end


import Base.==
function ==(ϕ₁::LatticeSite, ϕ₂::LatticeSite)
    return (ϕ₁.A₁ == ϕ₂.A₁ && ϕ₁.A₂ == ϕ₂.A₂ && ϕ₁.A₃ == ϕ₂.A₃ && ϕ₁.θ⁺ == ϕ₂.θ⁺ && 
        ϕ₁.θ⁻ == ϕ₂.θ⁻ && ϕ₁.u⁺ == ϕ₂.u⁺ && ϕ₁.u⁻ == ϕ₂.u⁻)
end

# -------------------------------------------------------------------------------------------------
# Utility function for running functions on single cuboids.
function getSubCuboidPropertyAux(funk::Function, chan::RemoteChannel{Channel{SubCuboid}})
    sc = fetch(chan)
    funk(sc)
end
function getSubCuboidProperty(funk::Function, chan::RemoteChannel{Channel{SubCuboid}})
    remotecall_fetch(getSubCuboidPropertyAux, chan.where, funk, chan)
end



############################################################################################################################
#                               Test functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

function testSCNeighbors(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    sc_pos = (1,1,1)
    display(sc.lattice[sc_pos[1]+1, sc_pos[2], sc_pos[3]] === sc.nb[sc_pos...].ϕᵣ₊₁)
    display(typeof(sc.shell))
    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        sc.nb[x,y,z].ϕᵣ₊₁.u⁺ = 1.0
    end
    put!(chan, sc)
    true
end

############################################################################################################################
#                               SubCuboid functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Ties a lattice, shell and constants together to form a sub-cuboid object by establishing the nearest neighbors of
# the lattice-sites on the lattice.
function SubCuboid(lattice::Array{LatticeSite, 3}, sc_nb::RemoteNeighbors{SubCuboid},
        shell::Vector{Array{LatticeSite, 2}}, shell_edge14::Vector{LatticeSite}, shell_edge16::Vector{LatticeSite},
        shell_edge23::Vector{LatticeSite}, shell_edge25::Vector{LatticeSite},
        shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite},
        consts::CubConstants, syst::SystConstants, sim::Controls, β::Float64)
    
    l₁ = consts.l₁; l₂ = consts.l₂; l₃ = consts.l₃
    (l₁ < 2 || l₂ < 2 || l₃ < 2) && throw(error("ERROR: Length < 2 when creating SubCuboid"))
    
    nb = latticeNeighbors(lattice, shell)
    nnb = latticeNextNeighbors(lattice, shell, shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45)
    SubCuboid(lattice, nb, nnb, sc_nb, shell, shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45,
              consts, syst, sim, β)
end

# Creates a SubCuboid and puts a reference to it in the remote-channel. This should be called from the process
# that the remote-channel is pointing to.
function remoteSubCuboid!(chan::RemoteChannel{Channel{SubCuboid}}, cub_lattice::Array{LatticeSite, 3}, 
        sc_nb::RemoteNeighbors{SubCuboid}, shell::Vector{Array{LatticeSite, 2}}, shell_edge14::Vector{LatticeSite},
        shell_edge16::Vector{LatticeSite}, shell_edge23::Vector{LatticeSite}, shell_edge25::Vector{LatticeSite},
        shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite},
        cub_consts::CubConstants, syst::SystConstants, sim::Controls, β::Float64)
    
    obj = SubCuboid(cub_lattice, sc_nb, shell, shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36,
                    shell_edge45, cub_consts, syst, sim, β)
    put!(chan, obj)
    nothing
end

# Given a grid of remote references to SubCuboids we construct the neighboring remote references to the position
# given by pos, given that we have periodic boundary conditions.
function RemoteNeighbors(grid::Array{RemoteChannel{Channel{SubCuboid}},3}, pos::Tuple{I,I,I}) where I<:Int
    s₁=size(grid,1); s₂=size(grid,2); s₃=size(grid,3)
    i₁,i₂,i₃ = pos
    nbᵣ₊₁ = grid[mod(i₁+1-1,s₁)+1, i₂, i₃]
    nbᵣ₋₁ = grid[mod(i₁-1-1,s₁)+1, i₂, i₃]
    nbᵣ₊₂ = grid[i₁, mod(i₂+1-1,s₂)+1, i₃]
    nbᵣ₋₂ = grid[i₁, mod(i₂-1-1,s₂)+1, i₃]
    nbᵣ₊₃ = grid[i₁, i₂, mod(i₃+1-1,s₃)+1]
    nbᵣ₋₃ = grid[i₁, i₂, mod(i₃-1-1,s₃)+1]
    nbᵣ₊₁₋₂ = grid[mod(i₁+1-1,s₁)+1, mod(i₂-1-1,s₂)+1, i₃]
    nbᵣ₊₁₋₃ = grid[mod(i₁+1-1,s₁)+1, i₂, mod(i₃-1-1,s₃)+1]
    nbᵣ₋₁₊₂ = grid[mod(i₁-1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃]
    nbᵣ₋₁₊₃ = grid[mod(i₁-1-1,s₁)+1, i₂, mod(i₃+1-1,s₃)+1]
    nbᵣ₊₂₋₃ = grid[i₁, mod(i₂+1-1,s₂)+1, mod(i₃-1-1,s₃)+1]
    nbᵣ₋₂₊₃ = grid[i₁, mod(i₂-1-1,s₂)+1, mod(i₃+1-1,s₃)+1]
    RemoteNeighbors{SubCuboid}(nbᵣ₊₁, nbᵣ₋₁, nbᵣ₊₂, nbᵣ₋₂, nbᵣ₊₃, nbᵣ₋₃, nbᵣ₊₁₋₂, nbᵣ₊₁₋₃, nbᵣ₋₁₊₂, nbᵣ₋₁₊₃, nbᵣ₊₂₋₃, nbᵣ₋₂₊₃)
end



############################################################################################################################
#                               Lattice Site Updating Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

#include("energy.jl")
include("xy_energy.jl")
include("functions_mc.jl")


############################################################################################################################
#                               Protocol Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

include("update_shell.jl")
include("update_protocol.jl")
