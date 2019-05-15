module CuboidModule

export SystConstants, Cuboid, mcSweep!, energy
export getShells, getShellEdges, shellSize, getLattice

using Distributed
using Distributions
using Test
using BenchmarkTools
#using Plots
#pyplot()

#@everywhere include("/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/sub_cuboid_mod.jl")
#src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
#@everywhere push!(LOAD_PATH, $src_path)
using SubCuboidModule



############################################################################################################################
#                               Types
#__________________________________________________________________________________________________________________________#
############################################################################################################################
#

# The State type keeps track of all the individual sub-cuboids that live in the various
# available processes
mutable struct Cuboid
    grid::Array{RemoteChannel{Channel{SubCuboid}}, 3}    # A 3D grid of the remote references to sub-cuboids.
    range_grid::Array{Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}, 3}
    syst::SystConstants
end



############################################################################################################################
#                               Construction Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Generate a lattice consisting of different choices for lattice sites
# 1:    Correlated mean field state with specified values for fields.
# 2:    Completely random state for all fields.
# 3:    u⁺ is random, the rest are mean field.
# 4:    u⁺ and u⁻ are random, the rest are mean field.
# 5:    θ⁺ and θ⁻ are random, the rest are mean field.
# 6:    All fields except u⁺ and u⁻ random.
function generateInitialLattice(choice::Int64, syst::SystConstants; u⁺=1.0, u⁻=0.0, θ⁺=0.0, θ⁻=0.0, A=[0.0, 0.0, 0.0])
    L₁ = syst.L₁; L₂ = syst.L₂; L₃ = syst.L₃
    (L₁ < 2 || L₂ < 2 || L₃ < 2) && throw(error("ERROR: Lattice size is too small for nearest neighbor interaction."))
    Amax = 10; umax = 3
    
    # Construct ordered state 
    if choice == 1
        # Mean field lattice sites for all fields.
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, u⁺, u⁻, x) for x=1:L₁, y=1:L₂, z=1:L₃]
    # Construct random state
    elseif choice == 2
        lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
                        rand(Uniform(0,2π)), rand(Uniform(0,2π)), umax*rand(), umax*rand(), x) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 3
        # Uniform mean field for all fields except u⁺ which is random.
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, umax*rand(), u⁻, x) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 4
        # Uniform mean field except for u⁺ and u⁻
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, umax*rand(), umax*rand(), x) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 5
        # Vary phases, all other fields are uniform mean fields.
        lattice = [LatticeSite([A[1], A[2], A[3]], rand(Uniform(0,2π)), rand(Uniform(0,2π)), u⁺, u⁻, x) 
            for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 6
        # All fields random except for amplitudes
        lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
                               rand(Uniform(0,2π)), rand(Uniform(0,2π)), u⁺, u⁻, x) for x=1:L₁, y=1:L₂, z=1:L₃]
    # We only have choices 1 - 6 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    lattice
end



# Given the lengths of dimensions of a 3D lattice and the number of splits of each dimension we create a
# corresponding 3D lattice of unit-range-triplets that consists of the ranges of each dimension of the
# original lattice that consistutes a sub-cube after the lattice is split.
function genSplittingRanges(Ls::Tuple{I,I,I}, splits::Tuple{I,I,I}) where I<:Int
    L₁, L₂, L₃ = Ls
    s₁, s₂, s₃ = splits
    
    # Make a grid of unit range triplets. These give the index extent of the sub-cuboid corresponding to the grid-point.
    l₁_min, n₁ = divrem(L₁,s₁)
    # n₁ is the number of grid-points having extra lattice points, while l₁_min is the minimum sub-cuboid length
    l₂_min, n₂ = divrem(L₂,s₂); l₃_min, n₃ = divrem(L₃,s₃);
    
    unit_range_grid = Array{Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}, 3}(undef, s₁,s₂,s₃)
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        if i₁ <= n₁
            range₁ = (i₁-1)*(l₁_min+1)+1:i₁*(l₁_min+1)
        else
            range₁ = n₁*(l₁_min+1) + (i₁-n₁-1)*l₁_min + 1:n₁*(l₁_min+1) + (i₁-n₁-1)*l₁_min + l₁_min
        end
        if i₂ <= n₂
            range₂ = (i₂-1)*(l₂_min+1)+1:i₂*(l₂_min+1)
        else
            range₂ = n₂*(l₂_min+1) + (i₂-n₂-1)*l₂_min + 1:n₂*(l₂_min+1) + (i₂-n₂-1)*l₂_min + l₂_min
        end
        if i₃ <= n₃
            range₃ = (i₃-1)*(l₃_min+1)+1:i₃*(l₃_min+1)
        else
            range₃ = n₃*(l₃_min+1) + (i₃-n₃-1)*l₃_min + 1:n₃*(l₃_min+1) + (i₃-n₃-1)*l₃_min + l₃_min
        end
        unit_range_grid[i₁, i₂, i₃] = (range₁, range₂, range₃)
    end
    unit_range_grid
end

# Given a splitting of the full lattice given by unit_range_grid, we construct a grid of shells that consists of the
# lattice sites at neighboring sub-cuboids.
function genShellGrid(lattice::Array{LatticeSite, 3}, s₁::Int64, s₂::Int64, s₃::Int64,
        unit_range_grid::Array{Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}, 3})
    
    shell_grid = Array{Vector{Array{LatticeSite, 2}}, 3}(undef, s₁,s₂,s₃)
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        # Get dimensions of sub-cuboid
        l₁ = length(unit_range_grid[i₁,i₂,i₃][1]); l₂ = length(unit_range_grid[i₁,i₂,i₃][2])
        l₃ = length(unit_range_grid[i₁,i₂,i₃][3])
        # Get the neighboring sub-cuboids (given periodic BC)
        cᵣ₊₁ = lattice[unit_range_grid[mod(i₁+1-1,s₁)+1, i₂, i₃]...]
        cᵣ₋₁ = lattice[unit_range_grid[mod(i₁-1-1,s₁)+1, i₂, i₃]...]
        cᵣ₊₂ = lattice[unit_range_grid[i₁, mod(i₂+1-1,s₂)+1, i₃]...]
        cᵣ₋₂ = lattice[unit_range_grid[i₁, mod(i₂-1-1,s₂)+1, i₃]...]
        cᵣ₊₃ = lattice[unit_range_grid[i₁, i₂, mod(i₃+1-1,s₃)+1]...]
        cᵣ₋₃ = lattice[unit_range_grid[i₁, i₂, mod(i₃-1-1,s₃)+1]...]
        shell = Array{Array{LatticeSite, 2}, 1}(undef, 6)
        shell[1] = cᵣ₊₁[1,:,:]
        shell[2] = cᵣ₋₁[end,:,:]
        shell[3] = cᵣ₊₂[:,1,:]
        shell[4] = cᵣ₋₂[:,end,:]
        shell[5] = cᵣ₊₃[:,:,1]
        shell[6] = cᵣ₋₃[:,:,end]
        # Save completed shell to shell_grid
        shell_grid[i₁,i₂,i₃] = shell
    end
    shell_grid
end

# Same as above but constructs one grid for each edge of the SubCuboid
function genShellEdgeGrid(lattice::Array{LatticeSite, 3}, s₁::I, s₂::I, s₃::I,
    range_grid::Array{Tuple{UnitRange{I},UnitRange{I},UnitRange{I}}, 3}) where I<:Int
    
    shell_edge14 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge16 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge23 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge25 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge36 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge45 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        # Get dimensions of sub-cuboid
        l₁ = length(range_grid[i₁,i₂,i₃][1]); l₂ = length(range_grid[i₁,i₂,i₃][2])
        l₃ = length(range_grid[i₁,i₂,i₃][3])
        # Get next-neighboring sub-cuboids (give periodic BC)
        scᵣ₊₁₋₂ = lattice[range_grid[mod(i₁+1-1,s₁)+1, mod(i₂-1-1,s₂)+1, i₃]...]
        scᵣ₊₁₋₃ = lattice[range_grid[mod(i₁+1-1,s₁)+1, i₂, mod(i₃-1-1,s₃)+1]...]
        scᵣ₋₁₊₂ = lattice[range_grid[mod(i₁-1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃]...]
        scᵣ₋₁₊₃ = lattice[range_grid[mod(i₁-1-1,s₁)+1, i₂, mod(i₃+1-1,s₃)+1]...]
        scᵣ₊₂₋₃ = lattice[range_grid[i₁, mod(i₂+1-1,s₂)+1, mod(i₃-1-1,s₃)+1]...]
        scᵣ₋₂₊₃ = lattice[range_grid[i₁, mod(i₂-1-1,s₂)+1, mod(i₃+1-1,s₃)+1]...]
        
        shell_edge14[i₁,i₂,i₃] = scᵣ₊₁₋₂[1,end,:]
        shell_edge16[i₁,i₂,i₃] = scᵣ₊₁₋₃[1,:,end]
        shell_edge23[i₁,i₂,i₃] = scᵣ₋₁₊₂[end,1,:]
        shell_edge25[i₁,i₂,i₃] = scᵣ₋₁₊₃[end,:,1]
        shell_edge36[i₁,i₂,i₃] = scᵣ₊₂₋₃[:,1,end]
        shell_edge45[i₁,i₂,i₃] = scᵣ₋₂₊₃[:,end,1]
    end
    
    shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45
end

# Given a full lattice with constants. Splits the lattice into sub-cuboids with corresponding nearest neighbor
# shells that live on separate processes.
function Cuboid(lattice::Array{LatticeSite, 3}, syst::SystConstants, sim::Controls, splits::Tuple{Int64, Int64, Int64}, β::Float64)
    s₁, s₂, s₃ = splits
    L₁ = size(lattice, 1); L₂ = size(lattice, 2); L₃ = size(lattice, 3)
    
    unit_range_grid = genSplittingRanges((L₁, L₂, L₃), splits)
    
    # Now we can make a sub-cuboid lattice with lattice[unit_range_grid[i₁, i₂, i₃]...]
    # We still need to make the shell around the sub-cuboid.
    
    shell_grid = genShellGrid(lattice, s₁,s₂,s₃, unit_range_grid)
    (shell_edge14_grid, shell_edge16_grid, shell_edge23_grid, shell_edge25_grid,
     shell_edge36_grid, shell_edge45_grid)  = genShellEdgeGrid(lattice, s₁,s₂,s₃, unit_range_grid)
    
    # With all the ingredients in hand, we just need to place the different sub-cuboids on the different available
    # worker processes.
    p = s₁*s₂*s₃ # The required number of processesRemoteChannel(()->Channel{Array{T,1}}(1), p)
    p <= nprocs()-1 || throw(error("ERROR: Not enough available workers to distribute sub-cubes"))
    grid = Array{RemoteChannel{Channel{SubCuboid}}, 3}(undef, s₁,s₂,s₃)
    
    # First we create remote channels in all grid points at pids starting on 2
    for i₃ = 1:s₃, i₂ = 1:s₂, i₁ = 1:s₁
        # First we create a remote-channel where we can put a sub-cube on a process
        pid = i₁ + (i₂-1)*s₁ + (i₃-1)*s₁*s₂ + 1 # We start on 2 and increase the pid
        grid[i₁,i₂,i₃] = RemoteChannel(()->Channel{SubCuboid}(1), pid)
    end
    
    for i₃ = 1:s₃, i₂ = 1:s₂, i₁ = 1:s₁
        pid = grid[i₁,i₂,i₃].where
        
        # Get neighboring remote channels
        sc_nb = RemoteNeighbors(grid, (i₁,i₂,i₃))
        
        # Then we create a SubCuboid and put it in the channel
        cub_lattice = lattice[unit_range_grid[i₁,i₂,i₃]...]
        shell = shell_grid[i₁,i₂,i₃]
        shell_edge14 = shell_edge14_grid[i₁,i₂,i₃]; shell_edge16 = shell_edge16_grid[i₁,i₂,i₃]
        shell_edge23 = shell_edge23_grid[i₁,i₂,i₃]; shell_edge25 = shell_edge25_grid[i₁,i₂,i₃]
        shell_edge36 = shell_edge36_grid[i₁,i₂,i₃]; shell_edge45 = shell_edge45_grid[i₁,i₂,i₃]
        l₁ = size(cub_lattice, 1); l₂ = size(cub_lattice, 2); l₃ = size(cub_lattice, 3)
        cub_consts = CubConstants(l₁, l₂, l₃)
        
        # We have the ingredients to make a SubCuboid object, but making it here would set the references to
        # be to the objects on this process, thus we need to create it on the worker-process
        
        remotecall_fetch(remoteSubCuboid!, pid, grid[i₁,i₂,i₃], cub_lattice, sc_nb, shell,
            shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45,
            cub_consts, syst, sim, β)
        # This calls the function SubCuboid on the process at pid, which should put a sub-cuboid in the
        # channel pointed to by the remote-channel grid[i₁,i₂,i₃] with correctly set neighbors.
        # A question here is if remotecall_fetch does a deep-copy of the function arguments such that
        # the shell lattice points will be copies of the original in the master process.
    end
    
    # Having distributed the sub-cuboid to their respective homes in the worker-processes, we collect their
    # references in a Cuboid object.
    Cuboid(grid, unit_range_grid, syst)
end

# Construct Cuboid using generation of lattice with options below.
# 1:    Correlated mean field state with specified values for fields.
# 2:    Completely random state for all fields.
# 3:    u⁺ is random, the rest are mean field.
# 4:    u⁺ and u⁻ are random, the rest are mean field.
# 5:    θ⁺ and θ⁻ are random, the rest are mean field.
# 6:    All fields except u⁺ and u⁻ random.
function Cuboid(choice::Int64, syst::SystConstants, splits::Tuple{I,I,I}, β::Float64;
                                sim=Controls(), u⁺=1.0, u⁻=0.0, θ⁺=0.0, θ⁻=0.0, A=[0.0, 0.0, 0.0]) where I<:Int
    lattice = generateInitialLattice(choice, syst; u⁺=u⁺, u⁻=u⁻, θ⁺=θ⁺, θ⁻=θ⁻, A=[A[1], A[2], A[3]])
    Cuboid(lattice, syst, sim, splits, β)
end


############################################################################################################################
#                               MC Functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# A function that goes through all LatticeSites in a Cuboid and follows the update-protocol
function mcSweep!(cub::Cuboid)
    # Create futures for all update-processes
    futs = [Future() for i = 1:length(cub.grid)]
    
    # Preform step 1 of updating internal points
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferInternalPoints!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 1")
    
    # Step 2: Updates intersection planes 4 and 1 in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionPlanes!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 2")
    
    # Step 3: Updates intersection lines 14 and neighboring two lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionLines!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 3")
    
    # Step 4: Updates plane 6 internal points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6Internals!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 4")
    
    # Step 5: Updates plane 6 intersection lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionLines!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 5")
    
    # Step 6: Updates intersection point and two neighboring plane 6 points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionPoint!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 6")

    nothing
end


############################################################################################################################
#                               Observables
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Calculates the energy for the whole cuboid in parallel
function energy(cub::Cuboid)
    E = 0.0
    futs = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where energy(chan)
    end
    for fut in futs
        E += fetch(fut)
    end
    E
end




############################################################################################################################
#                               Utility functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################
# Warning: slow function.
function getShells(rn::RemoteNeighbors{SubCuboid})
    scᵣ₊₁ = fetch(rn.rnᵣ₊₁)
    scᵣ₋₁ = fetch(rn.rnᵣ₋₁)
    scᵣ₊₂ = fetch(rn.rnᵣ₊₂)
    scᵣ₋₂ = fetch(rn.rnᵣ₋₂)
    scᵣ₊₃ = fetch(rn.rnᵣ₊₃)
    scᵣ₋₃ = fetch(rn.rnᵣ₋₃)
    
    shells = Array{Array{LatticeSite, 2},1}(undef, 6)
    shells[1] = scᵣ₊₁.lattice[1,:,:]
    shells[2] = scᵣ₋₁.lattice[end,:,:]
    shells[3] = scᵣ₊₂.lattice[:,1,:]
    shells[4] = scᵣ₋₂.lattice[:,end,:]
    shells[5] = scᵣ₊₃.lattice[:,:,1]
    shells[6] = scᵣ₋₃.lattice[:,:,end]
    shells
end

# Warning: slow function.
function getShellEdges(rn::SubCuboidModule.RemoteNeighbors{SubCuboid})
    scᵣ₊₁₋₂ = fetch(rn.rnᵣ₊₁₋₂)
    scᵣ₊₁₋₃ = fetch(rn.rnᵣ₊₁₋₃)
    scᵣ₋₁₊₂ = fetch(rn.rnᵣ₋₁₊₂)
    scᵣ₋₁₊₃ = fetch(rn.rnᵣ₋₁₊₃)
    scᵣ₊₂₋₃ = fetch(rn.rnᵣ₊₂₋₃)
    scᵣ₋₂₊₃ = fetch(rn.rnᵣ₋₂₊₃)
    
    edge14 = scᵣ₊₁₋₂.lattice[1, scᵣ₊₁₋₂.consts.l₂, :]
    edge16 = scᵣ₊₁₋₃.lattice[1, :, scᵣ₊₁₋₃.consts.l₃]
    edge23 = scᵣ₋₁₊₂.lattice[scᵣ₋₁₊₂.consts.l₁, 1, :]
    edge25 = scᵣ₋₁₊₃.lattice[scᵣ₋₁₊₃.consts.l₁, :, 1]
    edge36 = scᵣ₊₂₋₃.lattice[:, 1, scᵣ₊₂₋₃.consts.l₃]
    edge45 = scᵣ₋₂₊₃.lattice[:, scᵣ₋₂₊₃.consts.l₂, 1]
    
    edge14, edge16, edge23, edge25, edge36, edge45
end

# Function that collects the lattice from difference sub-cuboids and stitches it together to for a single lattice
# on the calling process. Note that this is a very slow process and should never be used in situations that are
# time-sensitive. This is mostly for debug purposes.
function getLattice(cub::Cuboid)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    
    lattice = Array{LatticeSite, 3}(undef, L₁, L₂, L₃)
    
    for (i, chan) = enumerate(cub.grid)
        sc = fetch(chan)
        lattice[cub.range_grid[i]...] .= sc.lattice
    end
    lattice
end

# Calculates the number of lattices sites that are stored in shells or shell_edges in total for all the
# Sub Cuboids in the Cuboid.
function shellSize(cub::Cuboid)
    si = 0
    for (i, ranges) = enumerate(cub.range_grid)
        l₁ = length(ranges[1]); l₂ = length(ranges[2]); l₃ = length(ranges[3])
        planes_size = 2*(l₁*l₂ + l₁*l₃ + l₂*l₃)
        edges_size = 2*(l₁ + l₂ + l₃)
        si += planes_size + edges_size
    end
    si
end

end
