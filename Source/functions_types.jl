#
#   Utility functions for the types
#
import Base.copy
function copy(ϕ::LatticeSite)
    LatticeSite([ϕ.A[1],ϕ.A[2]],ϕ.θ⁺,ϕ.θ⁻,ϕ.u⁺,ϕ.u⁻)
end
function copy(c::SystConstants)
    SystConstants(c.L, c.γ, c.g⁻², c.ν, c.f, c.β)
end
# Copy functions for Neighbors, NextNeighbors and NNNeighbors
function copy(nb::Neighbors)
	Neighbors(nb.ϕᵣ₊₁,nb.ϕᵣ₋₁,nb.ϕᵣ₊₂,nb.ϕᵣ₋₂)
end
function copy(nnb::NextNeighbors)
	NextNeighbors(nnb.ϕᵣ₊₁₊₂, nnb.ϕᵣ₊₁₋₂, nnb.ϕᵣ₋₁₊₂, nnb.ϕᵣ₋₁₋₂)
end
function copy(nnnb::NNNeighbors)
	NNNeighbors(nnnb.ϕᵣ₊₁₁, nnnb.ϕᵣ₋₁₁, nnnb.ϕᵣ₊₂₂, nnnb.ϕᵣ₋₂₂)
end
# Copy functions for LxL lattices of Neighbors
function copy(nbl::Array{Neighbors,2})
	L = size(nbl,2)
	nbl_copy = Array{Neighbors,2}(L,L)
	for h_pos = 1:L, v_pos = 1:L
		nbl_copy[v_pos,h_pos] = copy(nbl[v_pos,h_pos])
	end
	nbl_copy
end
function copy(nnbl::Array{NextNeighbors,2})
	L = size(nnbl,2)
	nnbl_copy = Array{NextNeighbors,2}(L,L)
	for h_pos = 1:L, v_pos = 1:L
		nnbl_copy[v_pos,h_pos] = copy(nnbl[v_pos,h_pos])
	end
	nnbl_copy
end
function copy(nnnbl::Array{NNNeighbors,2})
	L = size(nnnbl,2)
	nnnbl_copy = Array{NNNeighbors,2}(L,L)
	for h_pos = 1:L, v_pos = 1:L
		nnnbl_copy[v_pos,h_pos] = copy(nnnbl[v_pos,h_pos])
	end
	nnnbl_copy
end
function copy(ψ::State)
    Lx = size(ψ.lattice,2)
    Ly = size(ψ.lattice,1)
    lattice = [LatticeSite([ψ.lattice[y,x].A[1],ψ.lattice[y,x].A[2]],ψ.lattice[y,x].θ⁺,ψ.lattice[y,x].θ⁻,
            ψ.lattice[y,x].u⁺, ψ.lattice[y,x].u⁻) for y = 1:Ly, x=1:Lx]
	nbl = latticeNeighbors(lattice, Lx)
	nnbl = latticeNextNeighbors(lattice, Lx)
	nnnbl = latticeNNNeighbors(lattice, Lx)
    consts = copy(ψ.consts)
	State(lattice, consts, nbl, nnbl, nnnbl)
end
function copy(sim::Controls)
    return Controls(sim.θmax, sim.umax, sim.Amax)
end

# -------------------------------------------------------------------------------------------------
# Comparison operator for custom types

import Base.==
function ==(ϕ₁::LatticeSite, ϕ₂::LatticeSite)
    return (ϕ₁.A[1] == ϕ₂.A[1] && ϕ₁.A[2] == ϕ₂.A[2] && ϕ₁.θ⁺ == ϕ₂.θ⁺ && 
        ϕ₁.θ⁻ == ϕ₂.θ⁻ && ϕ₁.u⁺ == ϕ₂.u⁺ && ϕ₁.u⁻ == ϕ₂.u⁻)
end
function ==(c₁::SystConstants, c₂::SystConstants)
    return (c₁.L == c₂.L && c₁.γ == c₂.γ && c₁.g⁻² == c₂.g⁻² && c₁.ν == c₂.ν && c₁.f == c₂.f && c₁.β == c₂.β)
end
function ==(ψ₁::State, ψ₂::State)
    check = ψ₁.consts == ψ₂.consts
    L = ψ₁.consts.L
    for h_pos=1:L, v_pos=1:L
        check = check && ψ₁.lattice[v_pos,h_pos] == ψ₂.lattice[v_pos,h_pos]
    end
    return check
end

# -------------------------------------------------------------------------------------------------
# LatticeSite outer constructor
# Initializes a lattice site that has radom values
function LatticeSite()
    A_max = 3.0
    u⁺ = rand()
    LatticeSite([rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max, A_max))], two_pi*rand(), two_pi*rand(), 
        u⁺, √(1-u⁺^2))
end

# -------------------------------------------------------------------------------------------------
# Mutates target to have the same values as src.
function set!(target::LatticeSite, src::LatticeSite)
    target.A[1] = src.A[1]
    target.A[2] = src.A[2]
    target.θ⁺ = src.θ⁺
    target.θ⁻ = src.θ⁻
    target.u⁺ = src.u⁺
    target.u⁻ = src.u⁻
    return 1
end

####################################################################################################
#                            Functions for ::State
#
####################################################################################################

#TODO: Write lattice constructors for Neighbors, NextNeighbors and NNNeighbors

# -------------------------------------------------------------------------------------------------
# Outer constructor
# Initializes a state ψ that either 1: has zero as value for the the fluctuating gauge potential link variables,
# the phase and the u⁺ component (which means the u⁻=1) at each lattice site, or 2: has random values for
# these variables.
function State(choice::Int64, L::Int64)
    L <= 4 && throw(DomainError())

	consts = SystConstants(L, 1.0, 1/0.3^2, 0.3, 2.0/L, 0.5)
    # Construct ordered state 
    if choice == 1
        
        # Construct LxL lattice of LxL LatticeSites
        lattice = [LatticeSite([0,0],0,0,0,1) for y=1:L, x=1:L]
		nbl = latticeNeighbors(lattice,L)
		nnbl = latticeNextNeighbors(lattice,L)
		nnnbl = latticeNNNeighbors(lattice,L)
        ψ = State(lattice, consts, nbl, nnbl, nnnbl)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
						   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), 1) for y=1:L, x=1:L]
        for y=1:L, x=1:L
            lattice[y,x].u⁻ = √(1-lattice[y,x].u⁺^2)
        end
		nb = latticeNeighbors(lattice,L)
		nnb = latticeNextNeighbors(lattice,L)
		nnnb = latticeNNNeighbors(lattice,L)
        ψ = State(lattice, consts, nb, nnb, nnnb)
        
    # We only have choices 1 and 2 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    ψ
end
# Same as above but with non-default parameter-inputs from SystConstants
function State(choice::Int64, consts::SystConstants)
    N = consts.L
    N <= 4 && throw(DomainError())
    # Construct ordered state 
    if choice == 1
        
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([0,0],0,0,0,1) for y=1:N, x=1:N]
		nb = latticeNeighbors(lattice,N)
		nnb = latticeNextNeighbors(lattice,N)
		nnnb = latticeNNNeighbors(lattice,N)
        ψ = State(lattice, consts, nb, nnb, nnnb)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
						   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), 1) for y=1:N, x=1:N]
        for y=1:N, x=1:N
            lattice[y,x].u⁻ = √(1-lattice[y,x].u⁺^2)
        end
		nb = latticeNeighbors(lattice,N)
		nnb = latticeNextNeighbors(lattice,N)
		nnnb = latticeNNNeighbors(lattice,N)
        ψ = State(lattice, consts, nb, nnb, nnnb)
        
    # We only have choices 1 and 2 so far so other values for choice will give an error.
    else
        throw(DomainError())
    end
    ψ
end

# -------------------------------------------------------------------------------------------------
# Checks if any of the lattice sites in the lattices of two states are equal and returns false in
# this case, returns true otherwise.
function latticeCompletelyDifferent(ψ₁::State, ψ₂::State)
    ψ₁.consts.L == ψ₂.consts.L || throw(DomainError())
    isCompletelyDifferent = true
    for x=1:ψ₂.consts.L, y=1:ψ₂.consts.L
        isCompletelyDifferent = isCompletelyDifferent && (ψ₁.lattice[y,x].A[1] != ψ₂.lattice[y,x].A[1]
            && ψ₁.lattice[y,x].A[2] != ψ₂.lattice[y,x].A[2] && ψ₁.lattice[y,x].θ⁺ != ψ₂.lattice[y,x].θ⁺ 
            && ψ₁.lattice[y,x].θ⁻ != ψ₂.lattice[y,x].θ⁻ && ψ₁.lattice[y,x].u⁺ != ψ₂.lattice[y,x].u⁺)
    end
    return isCompletelyDifferent
end

# -------------------------------------------------------------------------------------------------
# Checks if the state is a physically possible state. Note that we have resitricted to only having
# LxL lattices.
function checkState(ψ::State)
    N = ψ.consts.L
    @test size(ψ.lattice, 1) == size(ψ.lattice, 2) == N
    @test N>2
    for x=1:N, y=1:N
        @test ψ.lattice[y,x].θ⁺ < 2π && ψ.lattice[y,x].θ⁺ >= 0
        @test ψ.lattice[y,x].θ⁻ < 2π && ψ.lattice[y,x].θ⁻ >= 0
        @test ψ.lattice[y,x].u⁺ <= 1 && ψ.lattice[y,x].u⁺ >= 0
        @test isapprox(ψ.lattice[y,x].u⁻^2+ψ.lattice[y,x].u⁺^2, 1.0, atol=0, rtol=1e-13)
    end
    @test typeof(ψ.consts.γ) == typeof(ψ.consts.g⁻²) == typeof(ψ.consts.ν) == typeof(ψ.consts.f) == Float64
	@test isapprox(ψ.consts.L*abs(ψ.consts.f) % 1, 0.0, atol=0, rtol=1e-13)
end

# -------------------------------------------------------------------------------------------------
# Saves the state to file with specified filename. If mode is set to "a" then the state is appended
# to the file.
function save(ψ::State, filename::AbstractString, mode::AbstractString="w")
    if mode != "w" && mode != "a"
        println("$(mode) is not a valid file-opening option.")
        return 0
    end
    println("Saving state to file $(filename)")
    L = ψ.consts.L
    open(filename, mode) do f
        # Write line indicating start of state
        write(f, "state start\n")
        # Writing system constants on top of the file
        write(f, "$(ψ.consts.L)\n")
        write(f, "$(ψ.consts.γ)\n")
        write(f, "$(ψ.consts.g⁻²)\n")
        write(f, "$(ψ.consts.ν)\n")
        write(f, "$(ψ.consts.f)\n")
        write(f, "$(ψ.consts.β)\n")
        # Then starting writing the state lattice
        for x=1:L, y=1:L
            ϕ = ψ.lattice[y,x]
            write(f, "$(ϕ.A[1]):$(ϕ.A[2]):$(ϕ.θ⁺):$(ϕ.θ⁻):$(ϕ.u⁺):$(ϕ.u⁻)\n")
        end
        write(f, "end\n")
    end
    return 1
end

# -------------------------------------------------------------------------------------------------
# The converse of the save function above: takes a text file where a state is saved and converts
# it into a State object which it return.
function readState(filename::AbstractString)
    open(filename, "r") do file
        line = readline(file)
        if line != "state start"
            throw(DomainError())
        end
        L = parse(Int64, readline(file))
        γ = parse(Float64, readline(file))
        g⁻² = parse(Float64, readline(file))
        ν = parse(Float64, readline(file))
        f = parse(Float64, readline(file))
        β = parse(Float64, readline(file))
        syst = SystConstants(L, γ, g⁻², ν, f, β)
        
        lattice = Array{LatticeSite, 2}(L,L)
        for h_pos = 1:L, v_pos = 1:L
            ϕ_values = split(readline(file), ":")
            A₁ = parse(Float64, ϕ_values[1])
            A₂ = parse(Float64, ϕ_values[2])
            θ⁺ = parse(Float64, ϕ_values[3])
            θ⁻ = parse(Float64, ϕ_values[4])
            u⁺ = parse(Float64, ϕ_values[5])
            u⁻ = parse(Float64, ϕ_values[6])
            lattice[v_pos,h_pos] = LatticeSite([A₁, A₂], θ⁺, θ⁻, u⁺, u⁻)
        end
        
        nbl = latticeNeighbors(lattice,L)
        nnbl = latticeNextNeighbors(lattice,L)
        nnnbl = latticeNNNeighbors(lattice,L)
        
        return State(lattice, syst, nbl, nnbl, nnnbl)
    end
end

####################################################################################################
#                            Functions for ::Controls
#
####################################################################################################

# -------------------------------------------------------------------------------------------------
function setValues!(sim_target::Controls, sim_source::Controls)
    sim_target.θmax = sim_source.θmax
    sim_target.umax = sim_source.umax
    sim_target.Amax = sim_source.Amax
    return
end


####################################################################################################
#                            Functions for ::Neighbors, ::NextNeighbors, ::NNNeighbors
#
####################################################################################################

# -------------------------------------------------------------------------------------------------