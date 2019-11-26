module CuboidModule

export SystConstants, Cuboid, mcSweep!, mcSweepEnUp!, energy, two_pi
export latticeMap, latticeSiteNeighborMap, latticeSiteMap
export shellSize, getLattice, fluxDensity, setTemp!, tuneUpdates!, printControls, estimateAR!
export getBeta, copy, avgZ, getTemp
export specificHeat, xyVortexSnapshot, vortexSnapshot, getSyst, getControls, chiralAmplitudeSnapshot
export retuneUpdates!, setUpdates!
export nMCS!, nMCSEnUp!, xyVortexSnapshotXYBasis
export twoCompHelMod, helMod, firstDerivativeTwist, secondDerivativeTwist

# For testing: compile with the exports below
export RemoteNeighbors, SubCuboid, LatticeSite

using Distributed
using Distributions
using Test
using BenchmarkTools
#using Plots
#pyplot()

#@everywhere include("/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/sub_cuboid_mod.jl")
#src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
#@everywhere push!(LOAD_PATH, $src_path)
#using SubCuboidModule
include("subcuboid_functions.jl")



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



###########################################################################################################################
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
        lattice = [LatticeSite(A[1], A[2], A[3], θ⁺, θ⁻, u⁺, u⁻, x, y) for x=1:L₁, y=1:L₂, z=1:L₃]
    # Construct random state
    elseif choice == 2
        lattice = [LatticeSite(rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),
                        rand(Uniform(0,2π))-π, rand(Uniform(0,2π))-π, umax*rand(), umax*rand(), x, y) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 3
        # Uniform mean field for all fields except u⁺ which is random.
        lattice = [LatticeSite(A[1], A[2], A[3], θ⁺, θ⁻, umax*rand(), u⁻, x, y) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 4
        # Uniform mean field except for u⁺ and u⁻
        lattice = [LatticeSite(A[1], A[2], A[3], θ⁺, θ⁻, umax*rand(), umax*rand(), x, y) for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 5
        # Vary phases, all other fields are uniform mean fields.
        lattice = [LatticeSite(A[1], A[2], A[3], rand(Uniform(0,2π))-π, rand(Uniform(0,2π))-π, u⁺, u⁻, x, y) 
            for x=1:L₁, y=1:L₂, z=1:L₃]
    elseif choice == 6
        # All fields random except for amplitudes
        lattice = [LatticeSite(rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),
                               rand(Uniform(0,2π))-π, rand(Uniform(0,2π))-π, u⁺, u⁻, x, y) for x=1:L₁, y=1:L₂, z=1:L₃]
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
   
    shell_edge13 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge14 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge16 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge23 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge24 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge25 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge36 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    shell_edge45 = Array{Vector{LatticeSite}, 3}(undef, s₁,s₂,s₃)
    
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        # Get dimensions of sub-cuboid
        l₁ = length(range_grid[i₁,i₂,i₃][1]); l₂ = length(range_grid[i₁,i₂,i₃][2])
        l₃ = length(range_grid[i₁,i₂,i₃][3])
        # Get next-neighboring sub-cuboids (give periodic BC)
        scᵣ₊₁₋₂ = lattice[range_grid[mod(i₁+1-1,s₁)+1, mod(i₂-1-1,s₂)+1, i₃]...]
        scᵣ₊₁₊₂ = lattice[range_grid[mod(i₁+1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃]...]
        scᵣ₊₁₋₃ = lattice[range_grid[mod(i₁+1-1,s₁)+1, i₂, mod(i₃-1-1,s₃)+1]...]
        scᵣ₋₁₊₂ = lattice[range_grid[mod(i₁-1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃]...]
        scᵣ₋₁₋₂ = lattice[range_grid[mod(i₁-1-1,s₁)+1, mod(i₂-1-1,s₂)+1, i₃]...]
        scᵣ₋₁₊₃ = lattice[range_grid[mod(i₁-1-1,s₁)+1, i₂, mod(i₃+1-1,s₃)+1]...]
        scᵣ₊₂₋₃ = lattice[range_grid[i₁, mod(i₂+1-1,s₂)+1, mod(i₃-1-1,s₃)+1]...]
        scᵣ₋₂₊₃ = lattice[range_grid[i₁, mod(i₂-1-1,s₂)+1, mod(i₃+1-1,s₃)+1]...]
       
        shell_edge13[i₁,i₂,i₃] = scᵣ₊₁₊₂[1,1,:]
        shell_edge14[i₁,i₂,i₃] = scᵣ₊₁₋₂[1,end,:]
        shell_edge16[i₁,i₂,i₃] = scᵣ₊₁₋₃[1,:,end]
        shell_edge23[i₁,i₂,i₃] = scᵣ₋₁₊₂[end,1,:]
        shell_edge24[i₁,i₂,i₃] = scᵣ₋₁₋₂[end,end,:]
        shell_edge25[i₁,i₂,i₃] = scᵣ₋₁₊₃[end,:,1]
        shell_edge36[i₁,i₂,i₃] = scᵣ₊₂₋₃[:,1,end]
        shell_edge45[i₁,i₂,i₃] = scᵣ₋₂₊₃[:,end,1]
    end
    
    shell_edge13, shell_edge14, shell_edge16, shell_edge23, shell_edge24, shell_edge25, shell_edge36, shell_edge45
end

# Given a full lattice with constants. Splits the lattice into sub-cuboids with corresponding nearest neighbor
# shells that live on separate processes.
function Cuboid(lattice::Array{LatticeSite, 3}, syst::SystConstants, sim::Controls, splits::Tuple{Int64, Int64, Int64}, β::Float64; pid_start=2)
    s₁, s₂, s₃ = splits
    L₁ = size(lattice, 1); L₂ = size(lattice, 2); L₃ = size(lattice, 3)
    
    unit_range_grid = genSplittingRanges((L₁, L₂, L₃), splits)
    
    # Now we can make a sub-cuboid lattice with lattice[unit_range_grid[i₁, i₂, i₃]...]
    # We still need to make the shell around the sub-cuboid.
    
    shell_grid = genShellGrid(lattice, s₁,s₂,s₃, unit_range_grid)
    (shell_edge13_grid, shell_edge14_grid, shell_edge16_grid, shell_edge23_grid, shell_edge24_grid, shell_edge25_grid,
     shell_edge36_grid, shell_edge45_grid)  = genShellEdgeGrid(lattice, s₁,s₂,s₃, unit_range_grid)
    
    # With all the ingredients in hand, we just need to place the different sub-cuboids on the different available
    # worker processes.
    p = s₁*s₂*s₃ # The required number of processesRemoteChannel(()->Channel{Array{T,1}}(1), p)
    pid_start-2 + p <= nprocs()-1 || throw(error("ERROR: Not enough available workers to distribute $(p) sub-cubes at pid $(pid_start) - $(pid_start + p-1)"))
    grid = Array{RemoteChannel{Channel{SubCuboid}}, 3}(undef, s₁,s₂,s₃)
    
    # First we create remote channels in all grid points at pids starting on pid_start
    for i₃ = 1:s₃, i₂ = 1:s₂, i₁ = 1:s₁
        # First we create a remote-channel where we can put a sub-cube on a process
        pid = i₁ + (i₂-1)*s₁ + (i₃-1)*s₁*s₂ + (pid_start-1) # We start on pid_start and increase the pid
        grid[i₁,i₂,i₃] = RemoteChannel(()->Channel{SubCuboid}(1), pid)
    end
    
    for i₃ = 1:s₃, i₂ = 1:s₂, i₁ = 1:s₁
        pid = grid[i₁,i₂,i₃].where
        
        # Get neighboring remote channels
        sc_nb = RemoteNeighbors(grid, (i₁,i₂,i₃))
        
        # Then we create a SubCuboid and put it in the channel
        cub_lattice = lattice[unit_range_grid[i₁,i₂,i₃]...]
        shell = shell_grid[i₁,i₂,i₃]
        shell_edge13 = shell_edge13_grid[i₁,i₂,i₃]; shell_edge14 = shell_edge14_grid[i₁,i₂,i₃]; shell_edge16 = shell_edge16_grid[i₁,i₂,i₃]
        shell_edge23 = shell_edge23_grid[i₁,i₂,i₃]; shell_edge24 = shell_edge24_grid[i₁,i₂,i₃]; shell_edge25 = shell_edge25_grid[i₁,i₂,i₃]
        shell_edge36 = shell_edge36_grid[i₁,i₂,i₃]; shell_edge45 = shell_edge45_grid[i₁,i₂,i₃]
        l₁ = size(cub_lattice, 1); l₂ = size(cub_lattice, 2); l₃ = size(cub_lattice, 3)
        cub_consts = CubConstants(l₁, l₂, l₃)
        
        # We have the ingredients to make a SubCuboid object, but making it here would set the references to
        # be to the objects on this process, thus we need to create it on the worker-process
        
        remotecall_fetch(remoteSubCuboid!, pid, grid[i₁,i₂,i₃], cub_lattice, sc_nb, shell,
            shell_edge13, shell_edge14, shell_edge16, shell_edge23, shell_edge24, shell_edge25,
            shell_edge36, shell_edge45, cub_consts, syst, sim, β)
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
                                sim=Controls(), u⁺=1.0, u⁻=0.0, θ⁺=0.0, θ⁻=0.0, A=[0.0, 0.0, 0.0], pid_start=2) where I<:Int
    lattice = generateInitialLattice(choice, syst; u⁺=u⁺, u⁻=u⁻, θ⁺=θ⁺, θ⁻=θ⁻, A=[A[1], A[2], A[3]])
    Cuboid(lattice, syst, sim, splits, β; pid_start=pid_start)
end

# Finds aviable process ids that are not occupied with sub-cuboids.
function findAvailProc(cubs::Array{Cuboid, 1})
    # Finding used workers
    used = Array{Int64, 1}(undef, 0)
    for cub in cubs
        for sc in cub.grid
            push!(used, sc.where)
        end
    end
    setdiff(vcat([1], workers()), used)
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
    
    # Step 2: Updates intersection plane 1 in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionPlanesStep2!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 2")

    # Step 3: Updates intersection plane 4 in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionPlanesStep3!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 3")

    # Step 4: Updates intersection lines 14 and neighboring two lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionLines!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 4")
    
    # Step 5: Updates plane 6 internal points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6Internals!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 5")
    
    # Step 6: Updates plane 6 intersection lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionLinesStep6!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 6")
    
    # Step 7: Updates plane 6 intersection lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionLinesStep7!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 7")
    
    # Step 8: Updates intersection point and two neighboring plane 6 points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionPoint!(chan)
    end
    for fut in futs; wait(fut); end
    #println("Step 8")

    nothing
end

function nMCS!(cub::Cuboid, n::Int)
    for i = 1:n; mcSweep!(cub); end
    nothing
end

# Runs n MCS on each of the states in the list in parallel. 
function nMCS!(cubs::Array{Cuboid, 1}, n::Int)
    futures = [Future() for cub in cubs]
    avail = findAvailProc(cubs)
#    start_pid = cubs[1].grid[1].where-length(cubs)
    for (i, cub) = enumerate(cubs)
        futures[i] = @spawnat cub.grid[end].where nMCS!(cub, n)
#        futures[i] = @spawnat avail[i] nMCS!(cub, n)
    end
    for fut in futures; wait(fut); end
    nothing
end


# Same as above but also return the energy difference and number of successful updates
function mcSweepEnUp!(cub::Cuboid)
    # Create futures for all update-processes
    futs = [Future() for i = 1:length(cub.grid)]
    en_diff = 0.0; updates = 0
    
    # Preform step 1 of updating internal points
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferInternalPointsRetEnUp!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end
    
    # Step 2: Updates intersection plane 1 in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionPlanesRetEnUpStep2!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end
    
    # Step 3: Updates intersection plan 4 in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionPlanesRetEnUpStep3!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end
    
    # Step 4: Updates intersection lines 14 and neighboring two lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferIntersectionLinesRetEnUp!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end
    
    # Step 5: Updates plane 6 internal points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6InternalsRetEnUp!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end
    
    # Step 6: Updates plane 6 intersection lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionLinesRetEnUpStep6!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end

    # Step 7: Updates plane 6 intersection lines in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionLinesRetEnUpStep7!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end

    
    # Step 8: Updates intersection point and two neighboring plane 6 points in parallel
    for (i, chan) = enumerate(cub.grid)
        futs[i] = @spawnat chan.where updateTransferPlane6IntersectionPointRetEnUp!(chan)
    end
    for fut in futs
        δE, up = fetch(fut)
        en_diff += δE; updates += up
    end

    en_diff, updates
end

function nMCSEnUp!(cub::Cuboid, n::Int)
    dEs = Array{Float64, 1}(undef, n)
    updates_list = Array{Int64, 1}(undef, n)
    for i = 1:n
        dEs[i], updates_list[i] = mcSweepEnUp!(cub)
    end
    dEs, updates_list
end

# Runs n mcSweepEnUp on each of the states in the list in parallel. 
function nMCSEnUp!(cubs::Array{Cuboid, 1}, n::Int)
    N_c = length(cubs)
    dE_lists = Array{Array{Float64, 1}, 1}(undef, N_c)
    updates_lists = Array{Array{Int64, 1}, 1}(undef, N_c)
    futures = [Future() for cub in cubs]
    avail = findAvailProc(cubs)
#    start_pid = cubs[1].grid[1].where-length(cubs)
    for (i, cub) = enumerate(cubs)
#        futures[i] = @spawnat (mod(start_pid+i-2,nprocs())+1) nMCSEnUp!(cub, n)
#        futures[i] = @spawnat avail[i] nMCSEnUp!(cub, n)
        futures[i] = @spawnat cub.grid[end].where nMCSEnUp!(cub, n)
    end
    for (i, fut) = enumerate(futures)
        dE_lists[i], updates_lists[i] = fetch(fut)
    end
    dE_lists, updates_lists
end

# -------------------------------------------------------------------------------------------------
# Takes a state cuboid and calculates estimates the acceptance rate using M
# Monte Carlo sweeps of the lattice. Then return the average and standard deviation of this
# fraction.
function estimateAR!(cub::Cuboid; M::Int64=500)
    fracs = zeros(M)
    N = cub.syst.L₁*cub.syst.L₂*cub.syst.L₃
    
    # Go through the entire lattice M times and gain the statistic of whether it gets updated or not
    for i = 1:M
        _, updated = mcSweepEnUp!(cub)
        fracs[i] = updated/N
    end
    return mean(fracs), std(fracs)
end

function estimateAR!(cubs::Array{Cuboid, 1}; M::Int64=500)
    Ns = [cub.syst.L₁*cub.syst.L₂*cub.syst.L₃ for cub in cubs]

    updated_lists = nMCSEnUp!(cubs, M)[2]
    fracs_list = Array{Array{Float64, 1}}(undef, length(cubs))
    for k = 1:length(cubs)
        fracs_list[k] = updated_lists[k]./Ns[k]
    end
    [(mean(fracs_list[k]), std(fracs_list[k])) for k = 1:length(cubs)]
end

############################################################################################################################
#                               Observables
#__________________________________________________________________________________________________________________________#
############################################################################################################################

include("helicity_modulus.jl")
include("XY_to_Chiral_transformation.jl")

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

# -----------------------------------------------------------------------------------------------------------
# Calculates the plaquette sum of the gauge field at a position pos, with the plaquette plane perpenducular
# to the x, y and z-axis. This corresponds to the different components of the curl of the gauge field
# in the continuum limit.
function fluxDensity(ϕ::LatticeSite, nb::CuboidModule.NearestNeighbors, s::SystConstants)
    ϕᵣ₊₁ = nb.ϕᵣ₊₁
    ϕᵣ₊₂ = nb.ϕᵣ₊₂
    ϕᵣ₊₃ = nb.ϕᵣ₊₃

    A₁, A₂, A₃ = linkVariables(ϕ, s)
    Aᵣ₊₂_₁ = linkVariableX(ϕᵣ₊₂, s)
    Aᵣ₊₂_₃ = linkVariableZ(ϕᵣ₊₂, s)
    Aᵣ₊₁_₂ = linkVariableY(ϕᵣ₊₁, s)
    Aᵣ₊₁_₃ = linkVariableZ(ϕᵣ₊₁, s)
    Aᵣ₊₃_₁ = linkVariableX(ϕᵣ₊₃, s)
    Aᵣ₊₃_₂ = linkVariableY(ϕᵣ₊₃, s)
    
    cur_A_x = (A₂ + Aᵣ₊₂_₃ - Aᵣ₊₃_₂ - A₃)
    cur_A_y = (A₃ + Aᵣ₊₃_₁ - Aᵣ₊₁_₃ - A₁)
    cur_A_z = (A₁ + Aᵣ₊₁_₂ - Aᵣ₊₂_₁ - A₂)
    
    cur_A_x, cur_A_y, cur_A_z
end
function fluxDensity(chan::RemoteChannel{Channel{SubCuboid}})
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    [fluxDensity(sc.lattice[x,y,z], sc.nb[x,y,z], sc.syst) for x = 1:l₁, y = 1:l₂, z = 1:l₃]
end
function fluxDensity(cub::Cuboid)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    
    ret_lattice = Array{Tuple{Float64,Float64,Float64},3}(undef, L₁, L₂, L₃)
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where fluxDensity(chan)
    end
    
    for (i, ranges) = enumerate(cub.range_grid)
        ret_lattice[ranges...] = fetch(futures[i])
    end
    
    ret_lattice
end


# Planar vorticity
# -----------------------------------------------------------------------------------------------------------
# New version of nᵣ based on suggestion by Troels. We are drawing the gauge-invariant phase
# difference back to [-π, π) instead of [0, 2π) and also adding 2πf so that we are measuring
# n instead of (n-f) which we would do by just doing the gauge-invariant phase difference.
function drawback(x::T) where T<:Real
    return mod2pi(x+π)-π
end
function nᵣ(syst::SystConstants, θᵣ₋₁::Float64, θᵣ::Float64, θᵣ₊₂::Float64, θᵣ₋₁₊₂::Float64, 
            A¹ᵣ₋₁::Float64, A²ᵣ::Float64, A¹ᵣ₋₁₊₂::Float64, A²ᵣ₋₁::Float64)
    (drawback(θᵣ - θᵣ₋₁ - A¹ᵣ₋₁) + drawback(θᵣ₊₂ - θᵣ - A²ᵣ)
        - drawback(θᵣ₊₂ - θᵣ₋₁₊₂ - A¹ᵣ₋₁₊₂)
        - drawback(θᵣ₋₁₊₂ - θᵣ₋₁  - A²ᵣ₋₁)+two_pi*syst.f)
end
#function nᵣ(c::SystConstants, ϕᵣ₋₁::LatticeSite, ϕᵣ::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, h_pos::Int64)
#    vort_θ⁺ = (drawback(ϕᵣ.θ⁺ - ϕᵣ₋₁.θ⁺ - ϕᵣ₋₁.A₁) + drawback(ϕᵣ₊₂.θ⁺ - ϕᵣ.θ⁺ - (ϕᵣ.A₂ + two_pi*c.f*(h_pos-1)))
#        - drawback(ϕᵣ₊₂.θ⁺ - ϕᵣ₋₁₊₂.θ⁺  - ϕᵣ₋₁₊₂.A₁)
#        - drawback(ϕᵣ₋₁₊₂.θ⁺ - ϕᵣ₋₁.θ⁺  - (ϕᵣ₋₁.A₂ + two_pi*c.f*(h_pos-2)))+two_pi*c.f)
#    vort_θ⁻ = (drawback(ϕᵣ.θ⁻ - ϕᵣ₋₁.θ⁻ - ϕᵣ₋₁.A₁) + drawback(ϕᵣ₊₂.θ⁻ - ϕᵣ.θ⁻ - (ϕᵣ.A₂ + two_pi*c.f*(h_pos-1)))
#        - drawback(ϕᵣ₊₂.θ⁻ - ϕᵣ₋₁₊₂.θ⁻  - ϕᵣ₋₁₊₂.A₁)
#        - drawback(ϕᵣ₋₁₊₂.θ⁻ - ϕᵣ₋₁.θ⁻  - (ϕᵣ₋₁.A₂ + two_pi*c.f*(h_pos-2)))+two_pi*c.f)
#    return vort_θ⁺, vort_θ⁻
#end
# Version of the above function which should be used when the energy used is in the xy-basis
function nᵣ(s::SystConstants, ϕᵣ₋₁::LatticeSite, ϕᵣ::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite)
    vort_θ⁺ = nᵣ(s, findθ⁺(ϕᵣ₋₁), findθ⁺(ϕᵣ), findθ⁺(ϕᵣ₊₂), findθ⁺(ϕᵣ₋₁₊₂), linkVariableX(ϕᵣ₋₁, s), linkVariableY(ϕᵣ, s),
                 linkVariableX(ϕᵣ₋₁₊₂, s), linkVariableY(ϕᵣ₋₁, s))
    vort_θ⁻ = nᵣ(s, findθ⁻(ϕᵣ₋₁), findθ⁻(ϕᵣ), findθ⁻(ϕᵣ₊₂), findθ⁻(ϕᵣ₋₁₊₂), linkVariableX(ϕᵣ₋₁, s), linkVariableY(ϕᵣ, s),
                 linkVariableX(ϕᵣ₋₁₊₂, s), linkVariableY(ϕᵣ₋₁, s))
    vort_θ⁺, vort_θ⁻
end
function nᵣ_XY(s::SystConstants, ϕᵣ₋₁::LatticeSite, ϕᵣ::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite)
    vort_θ⁺ = nᵣ(s, ϕᵣ₋₁.θ⁺, ϕᵣ.θ⁺, ϕᵣ₊₂.θ⁺, ϕᵣ₋₁₊₂.θ⁺, linkVariableX(ϕᵣ₋₁, s), linkVariableY(ϕᵣ, s),
                 linkVariableX(ϕᵣ₋₁₊₂, s), linkVariableY(ϕᵣ₋₁, s))
    vort_θ⁻ = nᵣ(s, ϕᵣ₋₁.θ⁻, ϕᵣ.θ⁻, ϕᵣ₊₂.θ⁻, ϕᵣ₋₁₊₂.θ⁻, linkVariableX(ϕᵣ₋₁, s), linkVariableY(ϕᵣ, s),
                 linkVariableX(ϕᵣ₋₁₊₂, s), linkVariableY(ϕᵣ₋₁, s))
    vort_θ⁺, vort_θ⁻
end

# -----------------------------------------------------------------------------------------------------------
# Given a lattice of lattice sites we find the corresponding vortex lattice assuming periodic BC.
# Warning: slow function. Should not be used to measure things for averages.
function vortexSnapshot(lattice::Array{LatticeSite, 3}, syst::SystConstants;
                        l₁ = size(lattice, 1), l₂ = size(lattice, 2), l₃ = size(lattice, 3), T = Float64)

    V⁺ = Array{T, 3}(undef, l₁,l₂,l₃); V⁻ = Array{T, 3}(undef, l₁,l₂,l₃)

    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ϕᵣ = lattice[x,y,z]
        ϕᵣ₋₁ = lattice[mod(x-1-1,l₁)+1,y,z]
        ϕᵣ₊₂ = lattice[x,mod(y+1-1,l₂)+1,z]
        ϕᵣ₋₁₊₂ = lattice[mod(x-1-1,l₁)+1,mod(y+1-1,l₂)+1,z]
        V⁺[x,y,z], V⁻[x,y,z] = nᵣ(syst, ϕᵣ₋₁, ϕᵣ, ϕᵣ₊₂, ϕᵣ₋₁₊₂,)
    end
    V⁺, V⁻
end

# --------------------------------------------------------------------------------------------------
# Projects 3D real space down on 2D through an averaging over the 3. dimension.
function avgZ(V::Array{T,3}) where T<:Real
    L₁ = size(V,1); L₂ = size(V, 2); L₃ = size(V,3)
    avg_V = zeros(T, L₁,L₂)
    for x = 1:L₁, y = 1:L₂
        for z = 1:L₃
            avg_V[x,y] += V[x,y,z]
        end
        avg_V[x,y] /= L₃
    end
    return avg_V
end
function sumZ(V::Array{T,3}) where T<:Real
    L₁ = size(V,1); L₂ = size(V, 2)
    reshape(sum(V; dims=3), L₁, L₂)
end

# -----------------------------------------------------------------------------------------------------------
# Assuming we have a state ψ, we want to find the lattice of vortices averaged over the z-direction.
# First we find the lattice of vorticies on a single SubCuboid then we return the average.
function xyVortexSnapshot(chan::RemoteChannel{Channel{SubCuboid}}; T=Float64)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    V⁺ = Array{T, 3}(undef, l₁,l₂,l₃); V⁻ = Array{T, 3}(undef, l₁,l₂,l₃)

    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ϕᵣ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        ϕᵣ₋₁ = nb.ϕᵣ₋₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂
        ϕᵣ₋₁₊₂ = sc.nnb[x,y,z].ϕᵣ₋₁₊₂
        V⁺[x,y,z], V⁻[x,y,z] = nᵣ(sc.syst, ϕᵣ₋₁, ϕᵣ, ϕᵣ₊₂, ϕᵣ₋₁₊₂)
    end

    sumZ(V⁺), sumZ(V⁻)
end

# Calculates a snapshot of the vortex lattice averaged over the z-direction. The vorticity is found in each
# lattice site with nᵣ∼ ∇×A, on each sub-cuboid, this lattice is summed over the z-direction, returned
# to the main-process which then averages these sums.
function xyVortexSnapshot(cub::Cuboid; T = Float64)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    s₁ = size(cub.grid, 1); s₂ = size(cub.grid, 2); s₃ = size(cub.grid, 3)
    
    V⁺ = Array{T, 3}(undef, L₁,L₂,s₃)
    V⁻ = Array{T, 3}(undef, L₁,L₂,s₃)

    # First we calculate the z-summed vortex lattice on each sub-cuboid.
    futures = [Future() for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃]
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        chan = cub.grid[i₁,i₂,i₃]
        futures[i₁,i₂,i₃] = @spawnat chan.where xyVortexSnapshot(chan)
    end
    
    # Each call to xyVortexSnapshot(::SubCuboid) returns an l₁×l₂ array whose position in the original lattice
    # is given by the range_grid
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        fut = futures[i₁,i₂,i₃]
        x_int, y_int, _ = cub.range_grid[i₁,i₂,i₃]

        V⁺[x_int,y_int,i₃], V⁻[x_int,y_int,i₃] = fetch(fut)
    end

    # Finally we calculate the average over the s₃ layers of averaged vorticity
    sumZ(V⁺)./L₃, sumZ(V⁻)./L₃
end

# Same two as above, but now for the native XY-basis.
function xyVortexSnapshotXYBasis(chan::RemoteChannel{Channel{SubCuboid}}; T=Float64)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    V⁺ = Array{T, 3}(undef, l₁,l₂,l₃); V⁻ = Array{T, 3}(undef, l₁,l₂,l₃)

    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ϕᵣ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        ϕᵣ₋₁ = nb.ϕᵣ₋₁
        ϕᵣ₊₂ = nb.ϕᵣ₊₂
        ϕᵣ₋₁₊₂ = sc.nnb[x,y,z].ϕᵣ₋₁₊₂
        V⁺[x,y,z], V⁻[x,y,z] = nᵣ_XY(sc.syst, ϕᵣ₋₁, ϕᵣ, ϕᵣ₊₂, ϕᵣ₋₁₊₂)
    end

    sumZ(V⁺), sumZ(V⁻)
end
function xyVortexSnapshotXYBasis(cub::Cuboid; T = Float64)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    s₁ = size(cub.grid, 1); s₂ = size(cub.grid, 2); s₃ = size(cub.grid, 3)
    
    V⁺ = Array{T, 3}(undef, L₁,L₂,s₃)
    V⁻ = Array{T, 3}(undef, L₁,L₂,s₃)

    # First we calculate the z-summed vortex lattice on each sub-cuboid.
    futures = [Future() for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃]
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        chan = cub.grid[i₁,i₂,i₃]
        futures[i₁,i₂,i₃] = @spawnat chan.where xyVortexSnapshotXYBasis(chan)
    end
    
    # Each call to xyVortexSnapshot(::SubCuboid) returns an l₁×l₂ array whose position in the original lattice
    # is given by the range_grid
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        fut = futures[i₁,i₂,i₃]
        x_int, y_int, _ = cub.range_grid[i₁,i₂,i₃]

        V⁺[x_int,y_int,i₃], V⁻[x_int,y_int,i₃] = fetch(fut)
    end

    # Finally we calculate the average over the s₃ layers of averaged vorticity
    sumZ(V⁺)./L₃, sumZ(V⁻)./L₃
end

# Amplitude measurements
# -----------------------------------------------------------------------------------------------------------

# Returns the chiral amplitudes of the system for each lattice site.
function chiralAmplitudeSnapshot(cub::Cuboid; T = Float64)
    u⁺_lattice = latticeSiteMap(findu⁺, cub, T)
    u⁻_lattice = latticeSiteMap(findu⁻, cub, T)

    u⁺_lattice, u⁻_lattice
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
function getShellEdges(rn::RemoteNeighbors{SubCuboid})
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

# The function is assumed to return the specified return type. The function is here
# assumed to take a SubCuboid object and position as arguments.
function latticeMap(funk::Function, chan::RemoteChannel{Channel{SubCuboid}}, T::DataType)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    ret_lattice = Array{T, 3}(undef, l₁,l₂,l₃)
    
    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ret_lattice[x,y,z] = T(funk(sc, (x,y,z)))
    end
    ret_lattice
end
# We assume that the function funk is defined on all processes. Call it on each sub-cuboid which in turn calls it on
# each lattice site, and stitch the results back into a lattice of function-returns.
function latticeMap(funk::Function, cub::Cuboid, T::DataType)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    
    ret_lattice = Array{T,3}(undef, L₁, L₂, L₃)
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where latticeMap(funk, chan, T)
    end
    
    for (i, ranges) = enumerate(cub.range_grid)
        ret_lattice[ranges...] = fetch(futures[i])
    end
    
    ret_lattice
end
# The function is here assumed to take a LatticeSite and a NearestNeighbor object as its two arguments and return
# the data-type specified by T
function latticeSiteNeighborMap(funk::Function, chan::RemoteChannel{Channel{SubCuboid}}, T::DataType)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    ret_lattice = Array{T, 3}(undef, l₁,l₂,l₃)
    
    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        ret_lattice[x,y,z] = T(funk(sc.lattice[x,y,z], sc.nb[x,y,z]))
    end
    ret_lattice
end
# We assume that the function funk is defined on all processes. Call it on each sub-cuboid which in turn calls it on
# each lattice site, and stitch the results back into a lattice of function-returns.
function latticeSiteNeighborMap(funk::Function, cub::Cuboid, T::DataType)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    
    ret_lattice = Array{T,3}(undef, L₁, L₂, L₃)
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where latticeSiteNeighborMap(funk, chan, T)
    end
    
    for (i, ranges) = enumerate(cub.range_grid)
        ret_lattice[ranges...] = fetch(futures[i])
    end
    
    ret_lattice
end
# The function is assumed only take a LatticeSite object as argument
function latticeSiteMap(funk::Function, chan::RemoteChannel{Channel{SubCuboid}}, T::DataType)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    [T(funk(sc.lattice[x,y,z])) for x = 1:l₁, y = 1:l₂, z = 1:l₃]
end
# We assume that the function funk is defined on all processes. Call it on each sub-cuboid which in turn calls it on
# each lattice site, and stitch the results back into a lattice of function-returns.
function latticeSiteMap(funk::Function, cub::Cuboid, T::DataType)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    
    ret_lattice = Array{T,3}(undef, L₁, L₂, L₃)
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where latticeSiteMap(funk, chan, T)
    end
    
    for (i, ranges) = enumerate(cub.range_grid)
        ret_lattice[ranges...] = fetch(futures[i])
    end
    
    ret_lattice
end

# The function funk is assumed to take a SubCuboid object as only argument.
function mutateSystem!(funk::Function, chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    funk(sc)
    put!(chan, sc)
    nothing
end
# We assume that the function funk is defined on all processes. The function should take a SubCuboid object as
# only argument and should manipulate it in some way.
function mutateSystem!(funk::Function, cub::Cuboid)
    
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where mutateSystem!(funk, chan)
    end
    
    for (i, ranges) = enumerate(cub.range_grid)
        wait(futures[i])
    end
    
    nothing
end

# The function is here assumed to take a LatticeSite and a NearestNeighbor object as its two first arguments
# as well as three ints coding for lattice position, and return
# the data-type specified by T. The sum of the values of the function is then returned.
function latticeSiteNeighborSum(funk::Function, chan::RemoteChannel{Channel{SubCuboid}}, ranges::Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}, T::DataType)
    sc = fetch(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    r₀₁ = ranges[1][1]-1; r₀₂ = ranges[2][1]-1; r₀₃ = ranges[3][1]-1

    s = T(0)
    
    for x = 1:l₁, y = 1:l₂, z = 1:l₃
        s += broadcast(T, funk(sc.lattice[x,y,z], sc.nb[x,y,z], r₀₁+x, r₀₂+y, r₀₃+z))
    end
    s
end
# We assume that the function funk is defined on all processes. Call it on each sub-cuboid which in turn calls it on
# each lattice site, then add the results of the function calls
function latticeSiteNeighborSum(funk::Function, cub::Cuboid, T::DataType)
    s = T(0)
    futures = [Future() for i = 1:length(cub.grid)]
    
    for (i, chan) = enumerate(cub.grid)
        futures[i] = @spawnat chan.where latticeSiteNeighborSum(funk, chan, cub.range_grid[i], T)
    end
    
    for i = 1:length(cub.grid)
        res = fetch(futures[i])
        s += T(res)
    end
    
    s
end


# Sets the temperature in each subcuboid of the cuboid.
function setTemp!(cub::Cuboid, T::R) where R<:Real
    mutateSystem!(sc -> sc.β = 1/T, cub)
    nothing
end

# Get the inverse temperature of the system
function getBeta(ψ::Cuboid)
    getSubCuboidProperty(sc -> sc.β, ψ.grid[1])
end
# Get the temperature
function getTemp(ψ::Cuboid)
    1/getBeta(ψ)
end

# Get the simulations constants in the system.
function getControls(cub::Cuboid)
    getSubCuboidProperty(sc -> sc.sim, cub.grid[1])
end

# Get the system constants in the system
function getSyst(cub::Cuboid)
    getSubCuboidProperty(sc -> sc.syst, cub.grid[1])
end

# -------------------------------------------------------------------------------------------------
# Sets the simulation constants of each subcuboid in the cuboid, either by specifying values
function setProposalIntervals!(cub::Cuboid; θ_max = π/3, u_max = 0.4, A_max = 3.0)
    mutateSystem!(sc -> sc.sim = Controls(θ_max, u_max, A_max), cub)
end
# or by specifying a ::Controls object
function setProposalIntervals!(cub::Cuboid, sim::Controls)
    mutateSystem!(sc -> sc.sim = sim, cub)
end

# -----------------------------------------------------------------------------------------------------------
# Adjusts the simulation controls such that it has an acceptance rate (AR) larger than LOWER.
# Tries to find the simulation controls that are simultaneously the ones with the highest values.
function tuneUpdates!(cub::Cuboid; M = 40, CUTOFF_MAX = 42, LOWER = 0.3,
        NEEDED_PROPOSALS = 5, DIVI = 1.5, TRIED_VALUES = 4)

    # M = 40                # How many MCS to do for estimating success fraction.
    # CUTOFF_MAX = 42       # How many times the while loop should run. 
    # LOWER = 0.3           # Minimum acceptance rate (AR). Must be less than 0.5
    # NEEDED_PROPOSALS = 5  # Number of sim with AR >= LOWER
    # DIVI = 1.5            # The number the current value is divided by to get lower limit of search interval.
    # TRIED_VALUES = 4      # Number of values to try in new interval
    
    sim = getControls(cub)
    proposedConstants = [copy(sim) for i=1:NEEDED_PROPOSALS]
    proposedAR = zeros(NEEDED_PROPOSALS)
    proposals = 0
    n = 0
    adjustment_mcs = 0
    
    # First we get an estimate of the accept probability
    (av, std) = estimateAR!(cub; M=M)
    adjustment_mcs = M
    if av >= LOWER
        return (av, adjustment_mcs, sim)
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
            setProposalIntervals!(cub, tries_A[i])
            (f_A[i], std) = estimateAR!(cub; M=M)
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
            setProposalIntervals!(cub, tries_u[i])
            (f_u[i], std) = estimateAR!(cub, M=M)
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
            setProposalIntervals!(cub, tries_θ[i])
            (f_θ[i], std) = estimateAR!(cub, M=M)
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
        i_max = argmax(accept_ratios)   # Finding index of sim that gave highest accept ratio.
        s₀ = tries[i_max]               # Setting this sim to the initial one.
        
        
        n += 1
        if n >= CUTOFF_MAX
            println("WARNING: Could not find simulation constant such that update probability 
                was higher than $(LOWER)")
            sim = s₀
            return (accept_ratios[i_max], adjustment_mcs, sim)
        end
    end
    # The end situation of the loop is that we have a number of proposals >= NEEDED_PROPOSALS.
    
    # Finding the distance from zero of the different proposals.
    norms = zeros(proposals)
    for i = 1:proposals
        norms[i] = proposedConstants[i].θmax^2 + proposedConstants[i].umax^2 + proposedConstants[i].Amax^2
    end
    i_max = argmax(norms)
    setValues!(sim, proposedConstants[i_max]) # Finally we update the simulation constants to 
    # the sim that has highest norm, and an accept ratio above LOWER.
    setProposalIntervals!(cub, sim)
    return (proposedAR[i_max], adjustment_mcs, sim)
    # Return the acceptance ratio of the new sim and the number of Monte-Carlo Sweeps done during this adjustment.
end

# -----------------------------------------------------------------------------------------------------------
# Tunes updates given manual starting parameters.
function retuneUpdates!(cub::Cuboid; θ_max = π/3, u_max = 0.4, A_max = 3.0)
    setProposalIntervals!(cub; θ_max = θ_max, u_max = u_max, A_max = A_max)
    tuneUpdates!(cub)
    nothing
end

# -----------------------------------------------------------------------------------------------------------
# Encapsulation for setProposalIntervals!
function setUpdates!(cub::Cuboid; θ_max = π/3, u_max = 0.4, A_max = 3.0)
    setProposalIntervals!(cub; θ_max = θ_max, u_max = u_max, A_max = A_max)
    nothing
end

# -------------------------------------------------------------------------------------------------
function printControls(cubs::Array{Cuboid, 1})
    println("State\tθmax\t\t\tumax\tAmax")
    for i = 1:length(cubs)
        sim = getControls(cubs[i])
        println("$i\t$(sim.θmax)\t$(sim.umax)\t$(sim.Amax)")
    end
    nothing
end
function printControls(sim::Controls)
    println("State\tθmax\t\t\tumax\tAmax")
    println("$i\t$(sim.θmax)\t$(sim.umax)\t$(sim.Amax)")
    nothing
end


import Base.copy
function copy(ψ::Cuboid)
    s₁ = size(ψ.grid, 1); s₂ = size(ψ.grid, 2); s₃ = size(ψ.grid, 3)
    split = (s₁, s₂, s₃)
    lattice = getLattice(ψ)
    sim = getControls(ψ)
    syst = getSyst(ψ)
    β = getBeta(ψ)
    pid_start = ψ.grid[1,1,1].where
    Cuboid(lattice, syst, sim, split, β; pid_start=pid_start)
end

end
