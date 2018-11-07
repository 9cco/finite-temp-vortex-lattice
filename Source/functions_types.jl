#
#   Utility functions for the types
#
import Base.copy
function copy(ϕ::LatticeSite)
    LatticeSite([ϕ.A[1],ϕ.A[2],ϕ.A[3]],ϕ.θ⁺,ϕ.θ⁻,ϕ.u⁺,ϕ.u⁻)
end
function copy(c::SystConstants)
    SystConstants(c.L, c.L₃, c.γ, c.g⁻², c.ν, c.κ₅, c.f, c.β)
end
# Copy functions for Neighbors, NextNeighbors and NNNeighbors
function copy(nb::NearestNeighbors)
	Neighbors(nb.ϕᵣ₊₁,nb.ϕᵣ₋₁,nb.ϕᵣ₊₂,nb.ϕᵣ₋₂)
end
function copy(nnb::NextNeighbors)
	NextNeighbors(nnb.ϕᵣ₊₁₊₂, nnb.ϕᵣ₊₁₋₂, nnb.ϕᵣ₋₁₊₂, nnb.ϕᵣ₋₁₋₂)
end
function copy(nnnb::NNNeighbors)
	NNNeighbors(nnnb.ϕᵣ₊₁₁, nnnb.ϕᵣ₋₁₁, nnnb.ϕᵣ₊₂₂, nnnb.ϕᵣ₋₂₂, nnnb.ϕᵣ₊₃₃, nnnb.ϕᵣ₋₃₃)
end
# Copy functions for LxL lattices of NearestNeighbors
function copy(nbl::Array{NearestNeighbors,2})
	L = size(nbl,2)
	nbl_copy = Array{NearestNeighbors,2}(L,L)
	for h_pos = 1:L, v_pos = 1:L
		nbl_copy[v_pos,h_pos] = copy(nbl[v_pos,h_pos])
	end
	nbl_copy
end
function copy(nbl::Array{NearestNeighbors,3})
	L = size(nbl,2)
    L₃ = size(nbl,3)
	nbl_copy = Array{NearestNeighbors,3}(L,L,L₃)
	for z_pos=1:L₃, h_pos = 1:L, v_pos = 1:L
		nbl_copy[v_pos,h_pos,z_pos] = copy(nbl[v_pos,h_pos,z_pos])
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
function copy(nnbl::Array{NextNeighbors,3})
	L = size(nnbl,2)
    L₃ = size(nnbl,3)
	nnbl_copy = Array{NextNeighbors,3}(L,L,L₃)
	for z_pos = 1:L₃, h_pos = 1:L, v_pos = 1:L
		nnbl_copy[v_pos,h_pos,z_pos] = copy(nnbl[v_pos,h_pos,z_pos])
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
function copy(nnnbl::Array{NNNeighbors,3})
	L = size(nnnbl,2)
    L₃ = size(nnnbl,3)
	nnnbl_copy = Array{NNNeighbors,3}(L,L,L₃)
	for z_pos = 1:L₃, h_pos = 1:L, v_pos = 1:L
		nnnbl_copy[v_pos,h_pos,z_pos] = copy(nnnbl[v_pos,h_pos,z_pos])
	end
	nnnbl_copy
end

function copy(ψ::State)
    Lx = size(ψ.lattice,2)
    Ly = size(ψ.lattice,1)
    L₃ = size(ψ.lattice,3)
    lattice = [copy(ψ.lattice[y,x,z]) for y=1:Ly, x=1:Lx, z=1:L₃]
	nbl = latticeNeighbors(lattice, Lx, L₃)
	nnbl = latticeNextNeighbors(lattice, Lx, L₃)
	nnnbl = latticeNNNeighbors(lattice, Lx, L₃)
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
    return (ϕ₁.A[1] == ϕ₂.A[1] && ϕ₁.A[2] == ϕ₂.A[2] && ϕ₁.A[3] == ϕ₂.A[3] && ϕ₁.θ⁺ == ϕ₂.θ⁺ && 
        ϕ₁.θ⁻ == ϕ₂.θ⁻ && ϕ₁.u⁺ == ϕ₂.u⁺ && ϕ₁.u⁻ == ϕ₂.u⁻)
end
function ==(c₁::SystConstants, c₂::SystConstants)
    return (c₁.L == c₂.L && c₁.L₃ == c₂.L₃ && c₁.γ == c₂.γ && c₁.g⁻² == c₂.g⁻² && c₁.ν == c₂.ν
            && c₁.κ₅ == c₂.κ₅ && c₁.f == c₂.f && c₁.β == c₂.β)
end
function ==(ψ₁::State, ψ₂::State)
    check = ψ₁.consts == ψ₂.consts
    for i in eachindex(ψ₁.lattice)
        check = check && ψ₁.lattice[i] == ψ₂.lattice[i]
    end
    return check
end

# -------------------------------------------------------------------------------------------------
# LatticeSite outer constructor
# Initializes a lattice site that has radom values
function LatticeSite()
    A_max = 3.0
    u⁺ = rand()
    LatticeSite([rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max, A_max)), rand(Uniform(-A_max,A_max))],
                two_pi*rand(), two_pi*rand(), u⁺, √(1-u⁺^2))
end

# -------------------------------------------------------------------------------------------------
# Mutates target to have the same values as src.
function set!(target::LatticeSite, src::LatticeSite)
    target.A[1] = src.A[1]
    target.A[2] = src.A[2]
    target.A[3] = src.A[3]
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

	consts = SystConstants(L, L, 1.0, 1/0.3^2, 0.3, 1.0, 2.0/L, 0.5)
    # Construct ordered state 
    if choice == 1
        
        # Construct LxL lattice of LxL LatticeSites
        lattice = [LatticeSite([0,0,0],0,0,0,1) for y=1:L, x=1:L, z=1:L]
		nbl = latticeNeighbors(lattice,L,L)
		nnbl = latticeNextNeighbors(lattice,L,L)
		nnnbl = latticeNNNeighbors(lattice,L,L)
        ψ = State(lattice, consts, nbl, nnbl, nnnbl)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
							  rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), rand()) for y=1:L, x=1:L, z=1:L]
#        for y=1:L, x=1:L
#            lattice[y,x].u⁻ = √(1-lattice[y,x].u⁺^2)
#        end
		nb = latticeNeighbors(lattice,L,L)
		nnb = latticeNextNeighbors(lattice,L,L)
		nnnb = latticeNNNeighbors(lattice,L,L)
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
    L₃ = consts.L₃
    N <= 4 && throw(DomainError())
    # Construct ordered state 
    if choice == 1
        
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([0,0,0],0,0,0,1) for y=1:N, x=1:N, z=1:L₃]
		nb = latticeNeighbors(lattice,N,L₃)
		nnb = latticeNextNeighbors(lattice,N,L₃)
		nnnb = latticeNNNeighbors(lattice,N,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^10
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
						   rand(Uniform(0,2π)), rand(Uniform(0,2π)), rand(), 1) for y=1:N, x=1:N, z=1:L₃]
        for y=1:N, x=1:N, z=1:L₃
            lattice[y,x,z].u⁻ = √(1-lattice[y,x,z].u⁺^2)
        end
		nb = latticeNeighbors(lattice,N,L₃)
		nnb = latticeNextNeighbors(lattice,N,L₃)
		nnnb = latticeNNNeighbors(lattice,N,L₃)
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
    ψ₁.consts.L₃ == ψ₂.consts.L₃ || throw(DomainError())
    isCompletelyDifferent = true
    for x=1:ψ₂.consts.L, y=1:ψ₂.consts.L, z = 1:ψ₂.consts.L₃
        isCompletelyDifferent = isCompletelyDifferent && (ψ₁.lattice[y,x,z].A[1] != ψ₂.lattice[y,x,z].A[1]
            && ψ₁.lattice[y,x,z].A[2] != ψ₂.lattice[y,x,z].A[2] && ψ₁.lattice[y,x,z].A[3] != ψ₂.lattice[y,x,z].A[3] 
            && ψ₁.lattice[y,x,z].θ⁺ != ψ₂.lattice[y,x,z].θ⁺ && ψ₁.lattice[y,x,z].θ⁻ != ψ₂.lattice[y,x,z].θ⁻ 
            && ψ₁.lattice[y,x,z].u⁺ != ψ₂.lattice[y,x,z].u⁺ && ψ₁.lattice[y,x,z].u⁻ != ψ₂.lattice[y,x,z].u⁻)
    end
    return isCompletelyDifferent
end

# -------------------------------------------------------------------------------------------------
# Checks if the state is a physically possible state. Note that we have resitricted to only having
# L×L×L₃ lattices.
function checkState(ψ::State)
    @test typeof(ψ.consts.L) == typeof(ψ.consts.L₃) == Int64
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    @test L > 0 && L₃ > 0
    @test size(ψ.lattice, 1) == size(ψ.lattice, 2) == L
    @test size(ψ.lattice, 3) == L₃
    @test L>2
    for x=1:L, y=1:L, z=1:L₃
        @test ψ.lattice[y,x,z].θ⁺ < 2π && ψ.lattice[y,x,z].θ⁺ >= 0
        @test ψ.lattice[y,x,z].θ⁻ < 2π && ψ.lattice[y,x,z].θ⁻ >= 0
        @test ψ.lattice[y,x,z].u⁺ >= 0
        @test ψ.lattice[y,x,z].u⁻ >= 0
#        @test isapprox(ψ.lattice[y,x].u⁻^2+ψ.lattice[y,x].u⁺^2, 1.0, atol=0, rtol=1e-13)
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
    L₃ = ψ.consts.L₃
    open(filename, mode) do f
        # Write line indicating start of state
        write(f, "state start\n")
        # Writing system constants on top of the file
        write(f, "$(ψ.consts.L)\n")
        write(f, "$(ψ.consts.L₃)\n")
        write(f, "$(ψ.consts.γ)\n")
        write(f, "$(ψ.consts.g⁻²)\n")
        write(f, "$(ψ.consts.ν)\n")
        write(f, "$(ψ.consts.κ₅)\n")
        write(f, "$(ψ.consts.f)\n")
        write(f, "$(ψ.consts.β)\n")
        # Then starting writing the state lattice
        for z=1:L₃, x=1:L, y=1:L
            ϕ = ψ.lattice[y,x,z]
            write(f, "$(ϕ.A[1]):$(ϕ.A[2]):$(ϕ.A[3]):$(ϕ.θ⁺):$(ϕ.θ⁻):$(ϕ.u⁺):$(ϕ.u⁻)\n")
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
        L₃ = parse(Int64, readline(file))
        γ = parse(Float64, readline(file))
        g⁻² = parse(Float64, readline(file))
        ν = parse(Float64, readline(file))
        κ₅ = parse(Float64, readline(file))
        f = parse(Float64, readline(file))
        β = parse(Float64, readline(file))
        syst = SystConstants(L, γ, g⁻², ν, f, β)
        
        lattice = Array{LatticeSite, 3}(L,L,L₃)
        for z_pos = 1:L₃, h_pos = 1:L, v_pos = 1:L
            ϕ_values = split(readline(file), ":")
            A₁ = parse(Float64, ϕ_values[1])
            A₂ = parse(Float64, ϕ_values[2])
            A₃ = parse(Float64, ϕ_values[3])
            θ⁺ = parse(Float64, ϕ_values[4])
            θ⁻ = parse(Float64, ϕ_values[5])
            u⁺ = parse(Float64, ϕ_values[6])
            u⁻ = parse(Float64, ϕ_values[7])
            lattice[v_pos,h_pos,z_pos] = LatticeSite([A₁, A₂, A₃], θ⁺, θ⁻, u⁺, u⁻)
        end
        
        nbl = latticeNeighbors(lattice,L)
        nnbl = latticeNextNeighbors(lattice,L)
        nnnbl = latticeNNNeighbors(lattice,L)
        
        return State(lattice, syst, nbl, nnbl, nnnbl)
    end
end

# -------------------------------------------------------------------------------------------------
# Return the average amplitudes in the system.
function meanAmplitudes(ψ::State)
	u⁺ = 0.0
	u⁻ = 0.0
    N = length(ψ.lattice)
    for ϕ in ψ.lattice
		u⁺ += ϕ.u⁺
		u⁻ += ϕ.u⁻
	end

	return u⁺/N, u⁻/N
end

# -------------------------------------------------------------------------------------------------
# Return the max and min amplitudes of the two components
# max(u⁺), min(u⁺), max(u⁻), min(u⁻)
function maxMinAmplitudes(ψ::State)
	ϕ = ψ.lattice[1,1]
	L = ψ.consts.L
	ex_u⁺ = [ϕ.u⁺, ϕ.u⁺]
	ex_u⁻ = [ϕ.u⁻, ϕ.u⁻]
    for ϕ in ψ.lattice
		if ϕ.u⁺ > ex_u⁺[1]
			ex_u⁺[1] = ϕ.u⁺
		elseif ϕ.u⁺ < ex_u⁺[2]
			ex_u⁺[2] = ϕ.u⁺
		end
		if ϕ.u⁻ > ex_u⁻[1]
			ex_u⁻[1] = ϕ.u⁻
		elseif ϕ.u⁻ < ex_u⁻[2]
			ex_u⁻[2] = ϕ.u⁻
		end
	end
	return ex_u⁺[1], ex_u⁺[2], ex_u⁻[1], ex_u⁻[2]
end

####################################################################################################
#                            Functions for ::Controls
#
####################################################################################################

function Controls()
    return Controls(π/3, 0.4, 3.0)
end

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
