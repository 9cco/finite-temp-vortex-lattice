module SubCuboidModule

using Distributed
using Distributions
using Test
using BenchmarkTools
#using Plots
#pyplot()

export remoteSubCuboid!, SubCuboid, RemoteNeighbors, NextNeighbors, NearestNeighbors, LatticeSite, SystConstants, Controls
export CubConstants, two_pi, copy, ==
export updateShell1!, updateShell2!, updateShell3!, updateShell4!, updateShell5!, updateShell6!, updateShellEdge14!, updateShellEdge23!
export internalPointsTransfer, intersectionPlanesTransfer, intersectionLinesTransfer, intersectionPointTransfer
export updateInternalPoints!, updateTransferInternalPoints!, updateIntersectionPlanes!, updateTransferIntersectionPlanes!
export updateIntersectionLines!, updateTransferIntersectionLines!, updateIntersectionPoint!, updateTransferIntersectionPoint!

#@everywhere struct Hack end
#function fixRC()
#    for p in workers()
#        @fetchfrom p Hack()
#    end
#end
#fixRC()

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
    A::Array{Float64,1}  # Fluctuating vector potential
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
	#ϕᵣ₊₁₊₂::LatticeSite
	ϕᵣ₊₁₋₂::LatticeSite
	ϕᵣ₋₁₊₂::LatticeSite
	#ϕᵣ₋₁₋₂::LatticeSite
    #ϕᵣ₊₁₋₃::LatticeSite
    #ϕᵣ₋₁₊₃::LatticeSite
    #ϕᵣ₊₂₋₃::LatticeSite
    #ϕᵣ₋₂₊₃::LatticeSite
end
struct RemoteNeighbors{T}
    rnᵣ₊₁::RemoteChannel{Channel{T}}
    rnᵣ₋₁::RemoteChannel{Channel{T}}
    rnᵣ₊₂::RemoteChannel{Channel{T}}
    rnᵣ₋₂::RemoteChannel{Channel{T}}
    rnᵣ₊₃::RemoteChannel{Channel{T}}
    rnᵣ₋₃::RemoteChannel{Channel{T}}
    rnᵣ₊₁₋₂::RemoteChannel{Channel{T}}
    rnᵣ₋₁₊₂::RemoteChannel{Channel{T}}
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

############################################################################################################################
#                               Neighbor functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Sets the nearest neighbors of the points in the lattice facing outwards from the sub-cube, excluding edges and corners,
# i.e. points that are internal to the periphery-planes of the sub-cube.
# A shell of nearest neighbor points is envisioned to encompass the lattice. This shell consists of 6 planes
# that are numbered so that 1: plane with x̂ normal, 2: plane with -x̂ normal, 3: plane with ŷ normal,
# 4: plane with -ŷ normal, 5: plane with ẑ normal, 6: plane with -ẑ normal.
function setBorderInternalPoints!(nb::Array{NearestNeighbors, 3}, lattice::Array{LatticeSite, 3}, 
        shell::Vector{Array{LatticeSite, 2}}, l₁::Int64, l₂::Int64, l₃::Int64)
    
    # Then we write neighbors for points in perifery planes that are internal to the planes
    for y = 2:l₂-1, z = 2:l₃-1
        # Plane 1 internal points
        ϕᵣ₊₁ = shell[1][y,z]
        ϕᵣ₋₁ = lattice[l₁-1, y, z]
        ϕᵣ₊₂ = lattice[l₁, y+1, z]
        ϕᵣ₋₂ = lattice[l₁, y-1, z]
        ϕᵣ₊₃ = lattice[l₁, y, z+1]
        ϕᵣ₋₃ = lattice[l₁, y, z-1]
        nb[l₁,y,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # Plane 2
        ϕᵣ₊₁ = lattice[2, y, z]
        ϕᵣ₋₁ = shell[2][y,z]
        ϕᵣ₊₂ = lattice[1, y+1, z]
        ϕᵣ₋₂ = lattice[1, y-1, z]
        ϕᵣ₊₃ = lattice[1, y, z+1]
        ϕᵣ₋₃ = lattice[1, y, z-1]
        nb[1,y,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    
    for x = 2:l₁-1, z = 2:l₃-1
        # Plane 3 internal points
        ϕᵣ₊₁ = lattice[x+1, l₂, z]
        ϕᵣ₋₁ = lattice[x-1, l₂, z]
        ϕᵣ₊₂ = shell[3][x,z]
        ϕᵣ₋₂ = lattice[x, l₂-1, z]
        ϕᵣ₊₃ = lattice[x, l₂, z+1]
        ϕᵣ₋₃ = lattice[x, l₂, z-1]
        nb[x,l₂,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # Plane 4 internal points
        ϕᵣ₊₁ = lattice[x+1, 1, z]
        ϕᵣ₋₁ = lattice[x-1, 1, z]
        ϕᵣ₊₂ = lattice[x, 2, z]
        ϕᵣ₋₂ = shell[4][x,z]
        ϕᵣ₊₃ = lattice[x, 1, z+1]
        ϕᵣ₋₃ = lattice[x, 1, z-1]
        nb[x,1,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    
    for x = 2:l₁-1, y = 2:l₂-1
        # Plane 5 internal points
        ϕᵣ₊₁ = lattice[x+1, y, l₃]
        ϕᵣ₋₁ = lattice[x-1, y, l₃]
        ϕᵣ₊₂ = lattice[x, y+1, l₃]
        ϕᵣ₋₂ = lattice[x, y-1, l₃]
        ϕᵣ₊₃ = shell[5][x,y]
        ϕᵣ₋₃ = lattice[x, y, l₃-1]
        nb[x,y,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # Plane 6 internal points
        ϕᵣ₊₁ = lattice[x+1, y, 1]
        ϕᵣ₋₁ = lattice[x-1, y, 1]
        ϕᵣ₊₂ = lattice[x, y+1, 1]
        ϕᵣ₋₂ = lattice[x, y-1, 1]
        ϕᵣ₊₃ = lattice[x, y, 2]
        ϕᵣ₋₃ = shell[6][x,y]
        nb[x,y,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    nothing
end

# Sets the nearest neighbors of the points in the lattice that are on the edges of the subcube (not including corners).
# A shell of nearest neighbor points is envisioned to encompass the lattice. This shell consists of 6 planes
# that are numbered so that 1: plane with x̂ normal, 2: plane with -x̂ normal, 3: plane with ŷ normal,
# 4: plane with -ŷ normal, 5: plane with ẑ normal, 6: plane with -ẑ normal.
function setBorderEdges!(nb::Array{NearestNeighbors, 3}, lattice::Array{LatticeSite, 3}, 
        shell::Vector{Array{LatticeSite, 2}}, l₁::Int64, l₂::Int64, l₃::Int64)
    
    # First we update vertical edges
    for z = 2:l₃-1
        # 41-edge (edge between planes 4 and 1), (x=l₁, y=1)
        ϕᵣ₊₁ = shell[1][1,z]
        ϕᵣ₋₁ = lattice[l₁-1, 1, z]
        ϕᵣ₊₂ = lattice[l₁, 2, z]
        ϕᵣ₋₂ = shell[4][l₁,z]
        ϕᵣ₊₃ = lattice[l₁, 1, z+1]
        ϕᵣ₋₃ = lattice[l₁, 1, z-1]
        nb[l₁,1,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 13-edge, (x=l₁, y=l₂)
        ϕᵣ₊₁ = shell[1][l₂,z]
        ϕᵣ₋₁ = lattice[l₁-1, l₂, z]
        ϕᵣ₊₂ = shell[3][l₁,z]
        ϕᵣ₋₂ = lattice[l₁, l₂-1, z]
        ϕᵣ₊₃ = lattice[l₁, l₂, z+1]
        ϕᵣ₋₃ = lattice[l₁, l₂, z-1]
        nb[l₁,l₂,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 32-edge (x=1, y=l₂)
        ϕᵣ₊₁ = lattice[2, l₂, z]
        ϕᵣ₋₁ = shell[2][l₂,z]
        ϕᵣ₊₂ = shell[3][1,z]
        ϕᵣ₋₂ = lattice[1, l₂-1, z]
        ϕᵣ₊₃ = lattice[1, l₂, z+1]
        ϕᵣ₋₃ = lattice[1, l₂, z-1]
        nb[1,l₂,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 24-edge (x=1, y=1)
        ϕᵣ₊₁ = lattice[2, 1, z]
        ϕᵣ₋₁ = shell[2][1,z]
        ϕᵣ₊₂ = lattice[1, 2, z]
        ϕᵣ₋₂ = shell[4][1,z]
        ϕᵣ₊₃ = lattice[1, 1, z+1]
        ϕᵣ₋₃ = lattice[1, 1, z-1]
        nb[1,1,z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    
    # Then we update the horizontal edges parallell to the x-axis
    for x = 2:l₁-1
        # 45-edge, (y=1, z=l₃)
        ϕᵣ₊₁ = lattice[x+1, 1, l₃]
        ϕᵣ₋₁ = lattice[x-1, 1, l₃]
        ϕᵣ₊₂ = lattice[x, 2, l₃]
        ϕᵣ₋₂ = shell[4][x,l₃]
        ϕᵣ₊₃ = shell[5][x,1]
        ϕᵣ₋₃ = lattice[x, 1, l₃-1]
        nb[x,1,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 53-edge, (y=l₂, z=l₃)
        ϕᵣ₊₁ = lattice[x+1, l₂, l₃]
        ϕᵣ₋₁ = lattice[x-1, l₂, l₃]
        ϕᵣ₊₂ = shell[3][x,l₃]
        ϕᵣ₋₂ = lattice[x, l₂-1, l₃]
        ϕᵣ₊₃ = shell[5][x,l₂]
        ϕᵣ₋₃ = lattice[x, l₂, l₃-1]
        nb[x,l₂,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 36-edge, (y=l₂, z=1)
        ϕᵣ₊₁ = lattice[x+1, l₂, 1]
        ϕᵣ₋₁ = lattice[x-1, l₂, 1]
        ϕᵣ₊₂ = shell[3][x,1]
        ϕᵣ₋₂ = lattice[x, l₂-1, 1]
        ϕᵣ₊₃ = lattice[x, l₂, 2]
        ϕᵣ₋₃ = shell[6][x,l₂]
        nb[x,l₂,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 64-edge, (y=1, z=1)
        ϕᵣ₊₁ = lattice[x+1, 1, 1]
        ϕᵣ₋₁ = lattice[x-1, 1, 1]
        ϕᵣ₊₂ = lattice[x, 2, 1]
        ϕᵣ₋₂ = shell[4][x,1]
        ϕᵣ₊₃ = lattice[x, 1, 2]
        ϕᵣ₋₃ = shell[6][x, 1]
        nb[x,1,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    
    # Then we update the horizontal edges parallell to the y-axis
    for y = 2:l₂-1
        # 51-edge, (x=l₁, z=l₃)
        ϕᵣ₊₁ = shell[1][y, l₃]
        ϕᵣ₋₁ = lattice[l₁-1, y, l₃]
        ϕᵣ₊₂ = lattice[l₁, y+1, l₃]
        ϕᵣ₋₂ = lattice[l₁, y-1, l₃]
        ϕᵣ₊₃ = shell[5][l₁, y]
        ϕᵣ₋₃ = lattice[l₁, y, l₃-1]
        nb[l₁,y,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 16-edge, (x=l₁, z=1)
        ϕᵣ₊₁ = shell[1][y, 1]
        ϕᵣ₋₁ = lattice[l₁-1, y, 1]
        ϕᵣ₊₂ = lattice[l₁, y+1, 1]
        ϕᵣ₋₂ = lattice[l₁, y-1, 1]
        ϕᵣ₊₃ = lattice[l₁, y, 2]
        ϕᵣ₋₃ = shell[6][l₁,y]
        nb[l₁,y,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 62-edge, (x=1, z=1)
        ϕᵣ₊₁ = lattice[2, y, 1]
        ϕᵣ₋₁ = shell[2][y, 1]
        ϕᵣ₊₂ = lattice[1, y+1, 1]
        ϕᵣ₋₂ = lattice[1, y-1, 1]
        ϕᵣ₊₃ = lattice[1, y, 2]
        ϕᵣ₋₃ = shell[6][1,y]
        nb[1,y,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
        
        # 25-edge, (x=1, z=l₃)
        ϕᵣ₊₁ = lattice[2, y, l₃]
        ϕᵣ₋₁ = shell[2][y, l₃]
        ϕᵣ₊₂ = lattice[1, y+1, l₃]
        ϕᵣ₋₂ = lattice[1, y-1, l₃]
        ϕᵣ₊₃ = shell[5][1,y]
        ϕᵣ₋₃ = lattice[1, y, l₃-1]
        nb[1,y,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    nothing
end

# Sets the nearest neighbors of the points in the lattice that are on the corners of the subcube.
# A shell of nearest neighbor points is envisioned to encompass the lattice. This shell consists of 6 planes
# that are numbered so that 1: plane with x̂ normal, 2: plane with -x̂ normal, 3: plane with ŷ normal,
# 4: plane with -ŷ normal, 5: plane with ẑ normal, 6: plane with -ẑ normal.
function setBorderCorners!(nb::Array{NearestNeighbors, 3}, lattice::Array{LatticeSite, 3}, 
        shell::Vector{Array{LatticeSite, 2}}, l₁::Int64, l₂::Int64, l₃::Int64)
    
    # Setting the 4 upper corners
    # 415-corner (corner intersected by the planes 4, 1 and 5), (x=l₁, y=1, z=l₃)
    ϕᵣ₊₁ = shell[1][1,l₃]
    ϕᵣ₋₁ = lattice[l₁-1, 1, l₃]
    ϕᵣ₊₂ = lattice[l₁, 2, l₃]
    ϕᵣ₋₂ = shell[4][l₁, l₃]
    ϕᵣ₊₃ = shell[5][l₁,1]
    ϕᵣ₋₃ = lattice[l₁, 1, l₃-1]
    nb[l₁,1,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 153-corner, (x=l₁, y=l₂, z=l₃)
    ϕᵣ₊₁ = shell[1][l₂,l₃]
    ϕᵣ₋₁ = lattice[l₁-1, l₂, l₃]
    ϕᵣ₊₂ = shell[3][l₁, l₃]
    ϕᵣ₋₂ = lattice[l₁, l₂-1, l₃]
    ϕᵣ₊₃ = shell[5][l₁,l₂]
    ϕᵣ₋₃ = lattice[l₁, l₂, l₃-1]
    nb[l₁,l₂,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 235-corner, (x=1, y=l₂, z=l₃)
    ϕᵣ₊₁ = lattice[2, l₂, l₃]
    ϕᵣ₋₁ = shell[2][l₂, l₃]
    ϕᵣ₊₂ = shell[3][1, l₃]
    ϕᵣ₋₂ = lattice[1, l₂-1, l₃]
    ϕᵣ₊₃ = shell[5][1,l₂]
    ϕᵣ₋₃ = lattice[1, l₂, l₃-1]
    nb[1,l₂,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 245-corner, (x=1, y=1, z=l₃)
    ϕᵣ₊₁ = lattice[2, 1, l₃]
    ϕᵣ₋₁ = shell[2][1, l₃]
    ϕᵣ₊₂ = lattice[1, 2, l₃]
    ϕᵣ₋₂ = shell[4][1,l₃]
    ϕᵣ₊₃ = shell[5][1,1]
    ϕᵣ₋₃ = lattice[1, 1, l₃-1]
    nb[1,1,l₃] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # Setting the 4 lower corners
    # 416-corner (corner intersected by the planes 4, 1 and 6), (x=l₁, y=1, z=1)
    ϕᵣ₊₁ = shell[1][1,1]
    ϕᵣ₋₁ = lattice[l₁-1, 1, 1]
    ϕᵣ₊₂ = lattice[l₁, 2, 1]
    ϕᵣ₋₂ = shell[4][l₁, 1]
    ϕᵣ₊₃ = lattice[l₁, 1, 2]
    ϕᵣ₋₃ = shell[6][l₁, 1]
    nb[l₁,1,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 163-corner, (x=l₁, y=l₂, z=1)
    ϕᵣ₊₁ = shell[1][l₂,1]
    ϕᵣ₋₁ = lattice[l₁-1, l₂, 1]
    ϕᵣ₊₂ = shell[3][l₁, 1]
    ϕᵣ₋₂ = lattice[l₁, l₂-1, 1]
    ϕᵣ₊₃ = lattice[l₁, l₂, 2]
    ϕᵣ₋₃ = shell[6][l₁, l₂]
    nb[l₁,l₂,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 236-corner, (x=1, y=l₂, z=1)
    ϕᵣ₊₁ = lattice[2, l₂, 1]
    ϕᵣ₋₁ = shell[2][l₂, 1]
    ϕᵣ₊₂ = shell[3][1, 1]
    ϕᵣ₋₂ = lattice[1, l₂-1, 1]
    ϕᵣ₊₃ = lattice[1, l₂, 2]
    ϕᵣ₋₃ = shell[6][1, l₂]
    nb[1,l₂,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    
    # 246-corner, (x=1, y=1, z=1)
    ϕᵣ₊₁ = lattice[2, 1, 1]
    ϕᵣ₋₁ = shell[2][1, 1]
    ϕᵣ₊₂ = lattice[1, 2, 1]
    ϕᵣ₋₂ = shell[4][1,1]
    ϕᵣ₊₃ = lattice[1, 1, 2]
    ϕᵣ₋₃ = shell[6][1, 1]
    nb[1,1,1] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    nothing
end

function latticeNeighbors(lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}})
    l₁ = size(lattice, 1); l₂ = size(lattice, 2); l₃ = size(lattice, 3)
    (l₁ < 1 || l₂ < 1 || l₃ < 1) && throw(DomainError())
    nb = Array{NearestNeighbors, 3}(undef, l₁, l₂, l₃)
    
    # First we write the trivial nearest neighbors of internal points
    for x = 2:l₁-1, y = 2:l₂-1, z = 2:l₃-1
        ϕᵣ₊₁ = lattice[x+1, y, z]
        ϕᵣ₋₁ = lattice[x-1, y, z]
        ϕᵣ₊₂ = lattice[x, y+1, z]
        ϕᵣ₋₂ = lattice[x, y-1, z]
        ϕᵣ₊₃ = lattice[x, y, z+1]
        ϕᵣ₋₃ = lattice[x, y, z-1]
        nb[x, y, z] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    
    # Set the nearest neighbors of points in the periphery of the lattice, except for edges and corners.
    setBorderInternalPoints!(nb, lattice, shell, l₁, l₂, l₃)
    
    # Set the nearest neighbors of edge-points in the lattice , except for corners.
    setBorderEdges!(nb, lattice, shell, l₁, l₂, l₃)
    
    # Set the nearest neighbors of corner-points on the lattice
    setBorderCorners!(nb, lattice, shell, l₁, l₂, l₃)
    nb
end

# ---------------------------------------------------------------------------------------------------
# Set the next nearest neighbors on the lattice of a sub-cuboid with shells around it.
function latticeNextNeighbors(lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}},
    shell_edge14::Vector{LatticeSite}, shell_edge23::Vector{LatticeSite})
    l₁ = size(lattice, 1); l₂ = size(lattice, 2); l₃ = size(lattice, 3)
    (l₁ < 1 || l₂ < 1 || l₃ < 1) && throw(DomainError())
    nnb = Array{NextNeighbors, 3}(undef, l₁, l₂, l₃)
    
    # First we do all elementary internal xy-plane points
    for x = 2:l₁-1, y = 2:l₂-1, z = 1:l₃
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
    end
    
    # Next we do the upper-edge xy-plane points, which are in border-plane 3. (y=l₂)
    for x = 2:l₁-1, z = 1:l₃
        ϕᵣ₊₁₋₂ = lattice[x+1,l₂-1,z]
        ϕᵣ₋₁₊₂ = shell[3][x-1,z]
        nnb[x,l₂,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
        
        # Next we do the lower-edge xy-plane points, where are in border-plane 4, (y=1)
        ϕᵣ₊₁₋₂ = shell[4][x+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,2,z]
        nnb[x,1,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
    end
    
    # Next we do the left-edge xy-plane points, which are in border-plane 2, (x=1)
    for y = 2:l₂-1, z = 1:l₃
        ϕᵣ₊₁₋₂ = lattice[2,y-1,z]
        ϕᵣ₋₁₊₂ = shell[2][y+1,z]
        nnb[1,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
        
        # Next we do the right-edge xy-plane points, which er in border-plane 1, (x=l₁)
        ϕᵣ₊₁₋₂ = shell[1][y-1,z]
        ϕᵣ₋₁₊₂ = lattice[l₁-1,y+1,z]
        nnb[l₁,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
    end
    
    # Finally we do the edges of the cube which are the corner xy-plane points.
    for z = 1:l₃
        # 23-edge (x=1,y=l₂)
        ϕᵣ₊₁₋₂ = lattice[2,l₂-1,z]
        ϕᵣ₋₁₊₂ = shell_edge23[z]
        nnb[1,l₂,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
        
        # 14-edge (x=l₁,y=1)
        ϕᵣ₊₁₋₂ = shell_edge14[z]
        ϕᵣ₋₁₊₂ = lattice[l₁-1,1+1,z]
        nnb[l₁,1,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
        
        # 24-edge (x=1,y=1)
        ϕᵣ₊₁₋₂ = shell[4][1+1,z]
        ϕᵣ₋₁₊₂ = shell[2][1+1,z]
        nnb[1,1,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
        
        # 13-edge (x=l₁,y=l₂)
        ϕᵣ₊₁₋₂ = shell[1][l₂-1,z]
        ϕᵣ₋₁₊₂ = shell[3][l₁-1,z]
        nnb[l₁,l₂,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂)
    end
    nnb
end


############################################################################################################################
#                               Utility functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

import Base.copy
function copy(re_nb::RemoteNeighbors{T}) where T
    rnᵣ₊₁ = re_nb.rnᵣ₊₁
    rnᵣ₋₁ = re_nb.rnᵣ₋₁
    rnᵣ₊₂ = re_nb.rnᵣ₊₂
    rnᵣ₋₂ = re_nb.rnᵣ₋₂
    rnᵣ₊₃ = re_nb.rnᵣ₊₃
    rnᵣ₋₃ = re_nb.rnᵣ₋₃
    rnᵣ₊₁₋₂ = re_nb.rnᵣ₊₁₋₂
    rnᵣ₋₁₊₂ = re_nb.rnᵣ₋₁₊₂
    RemoteNeighbors{T}(rnᵣ₊₁, rnᵣ₋₁, rnᵣ₊₂, rnᵣ₋₂, rnᵣ₊₃, rnᵣ₋₃, rnᵣ₊₁₋₂, rnᵣ₋₁₊₂)
end

import Base.==
function ==(ϕ₁::LatticeSite, ϕ₂::LatticeSite)
    return (ϕ₁.A[1] == ϕ₂.A[1] && ϕ₁.A[2] == ϕ₂.A[2] && ϕ₁.A[3] == ϕ₂.A[3] && ϕ₁.θ⁺ == ϕ₂.θ⁺ && 
        ϕ₁.θ⁻ == ϕ₂.θ⁻ && ϕ₁.u⁺ == ϕ₂.u⁺ && ϕ₁.u⁻ == ϕ₂.u⁻)
end

############################################################################################################################
#                               SubCuboid functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Ties a lattice, shell and constants together to form a sub-cuboid object by establishing the nearest neighbors of
# the lattice-sites on the lattice.
function SubCuboid(lattice::Array{LatticeSite, 3}, sc_nb::RemoteNeighbors{SubCuboid},
        shell::Vector{Array{LatticeSite, 2}}, shell_edge14::Vector{LatticeSite}, shell_edge23::Vector{LatticeSite}, 
        consts::CubConstants, syst::SystConstants, sim::Controls, β::Float64)
    
    l₁ = consts.l₁; l₂ = consts.l₂; l₃ = consts.l₃
    (l₁ < 2 || l₂ < 2 || l₃ < 2) && throw(error("ERROR: Length < 2 when creating SubCuboid"))
    
    nb = latticeNeighbors(lattice, shell)
    nnb = latticeNextNeighbors(lattice, shell, shell_edge14, shell_edge23)
    SubCuboid(lattice, nb, nnb, sc_nb, shell, shell_edge14, shell_edge23, consts, syst, sim, β)
end

# Creates a SubCuboid and puts a reference to it in the remote-channel. This should be called from the process
# that the remote-channel is pointing to.
function remoteSubCuboid!(chan::RemoteChannel{Channel{SubCuboid}}, cub_lattice::Array{LatticeSite, 3}, 
        sc_nb::RemoteNeighbors{SubCuboid}, shell::Vector{Array{LatticeSite, 2}}, shell_edge14::Vector{LatticeSite}, 
        shell_edge23::Vector{LatticeSite}, cub_consts::CubConstants, syst::SystConstants, sim::Controls, β::Float64)
    
    obj = SubCuboid(cub_lattice, sc_nb, shell, shell_edge14, shell_edge23, cub_consts, syst, sim, β)
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
    nbᵣ₋₁₊₂ = grid[mod(i₁-1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃]
    RemoteNeighbors{SubCuboid}(nbᵣ₊₁, nbᵣ₋₁, nbᵣ₊₂, nbᵣ₋₂, nbᵣ₊₃, nbᵣ₋₃, nbᵣ₊₁₋₂, nbᵣ₋₁₊₂)
end


############################################################################################################################
#                               Update Shell functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# These types of functions are typically initiated through a remote-call from a different process containing a sub-cuboid
# whose lattice has just been updated and thus the function calls are ment to reflect this change in the neighboring
# sub-cuboid.

# Takes the remote channel of a sub-cuboid and updates its shell 1 with updated values in shell_update.
# shell_update is encoded with the y,z-position in the shell, followed by the updated value.
function updateShell1!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell1 = sc.shell[1]
    for update in shell_update
        y, z, ϕ = update
        shell1[y,z] = ϕ
    end
    put!(chan, sc)
    nothing
end

# Takes the remote channel of a sub-cuboid and updates its shell 2 with updated values in shell_update.
# shell_update is encoded with the y,z-position in the shell, followed by the updated value.
function updateShell2!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell2 = sc.shell[2]
    for update in shell_update
        y, z, ϕ = update
        shell2[y,z] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShell3!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell3 = sc.shell[3]
    for update in shell_update
        x, z, ϕ = update
        shell3[x,z] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShell4!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell4 = sc.shell[4]
    for update in shell_update
        x, z, ϕ = update
        shell4[x,z] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShell5!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell5 = sc.shell[5]
    for update in shell_update
        x, y, ϕ = update
        shell5[x,y] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShell6!(chan::RemoteChannel{Channel{SubCuboid}}, 
        shell_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    shell6 = sc.shell[6]
    for update in shell_update
        x, y, ϕ = update
        shell6[x,y] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShellEdge14!(chan::RemoteChannel{Channel{SubCuboid}}, 
        edge_update::Vector{Tuple{I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    for update in edge_update
        z, ϕ = update
        sc.shell_edge14[z] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShellEdge16!(chan::RemoteChannel{Channel{SubCuboid}},
                            edge_update::Vector{Tuple{I,LatticeSite}}) where I<:Int

    sc = take!(chan)
    for update in edge_update
        y, ϕ = update
        sc.shell_edge16[z] = ϕ
    end
    put!(chan, sc)
    nothing
end
function updateShellEdge23!(chan::RemoteChannel{Channel{SubCuboid}}, 
        edge_update::Vector{Tuple{I,LatticeSite}}) where I<:Int
    
    sc = take!(chan)
    for update in edge_update
        z, ϕ = update
        sc.shell_edge23[z] = ϕ
    end
    put!(chan, sc)
    nothing
end



############################################################################################################################
#                               Transfer functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################
#
# Then we move up one abstraction layer and use these update functions in functions that transfers information from
# a sub-cuboid that has had certain parts of it updated. The functions below is ment to be used to implement the
# nearest-neighbor parallel sub-cuboid update protocol, which makes sure that no two lattice points that have a nearest
# neighbor interaction or ϕᵣ₋₁₊₂, ϕᵣ₊₁₋₂ interaction and are on different sub-cuboids are updated simultaneously,
# but instead are updated one after the other.


# Suppose we have updated the 3 points (x=l₁-1, y=1), (x=l₁, y=1) and (x=l₁, y=2) in the z=1 plane. This function sends
# information to the shells and shell-edges of neighboring sub-cuboids to reflect this.
function intersectionPointTransfer(sc_nb::RemoteNeighbors{SubCuboid}, 
        plane1_update::Vector{Tuple{I,I,LatticeSite}}, plane4_update::Vector{Tuple{I,I,LatticeSite}},
        plane6_update::Vector{Tuple{I,I,LatticeSite}}, edge14_update::Vector{Tuple{I,LatticeSite}}) where I<:Int
    
    # SubCuboid scᵣ₋₂ is affected by plane 4 updates through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, plane4_update)
    
    # SubCuboid sᵣ₊₁ is affected by plane 1 updates through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, plane1_update)

    # SubCuboid scᵣ₋₃ is affected by plane 6 updates through its shell[5]
    chan = sc_nb.rnᵣ₋₃
    remotecall_wait(updateShell5!, chan.where, chan, plane6_update)
    
    # SubCuboid scᵣ₊₁₋₂ is affected by the intersection point update in edge14_update through its shell_edge23
    chan = sc_nb.rnᵣ₊₁₋₂
    remotecall_wait(updateShellEdge23!, chan.where, chan, edge14_update)
    nothing
end

############################################################################################################################
#                               Lattice Site Updating Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

include("functions_mc.jl")


############################################################################################################################
#                               Protocol Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

include("functions_update_protocol.jl")


end
