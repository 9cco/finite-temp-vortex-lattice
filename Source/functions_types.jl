#
#   Utility functions for the types
#

####################################################################################################
#                            Functions for ::Controls
#
####################################################################################################

function Controls(θ_max::Float64, u_max::Float64, A_max::Float64)
    return Controls(θ_max, u_max, A_max, Uniform(-θ_max,θ_max), Uniform(-u_max,u_max), Uniform(-A_max, A_max))
end

function Controls()
    return Controls(π/3, 0.4, 3.0)
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

# -------------------------------------------------------------------------------------------------
function printSimControls(sim_list::Array{Controls})
    println("State\tθmax\t\t\tumax\tAmax")
    for i = 1:length(sim_list)
        println("$i\t$(sim_list[i].θmax)\t$(sim_list[i].umax)\t$(sim_list[i].Amax)")
    end
    return
end

####################################################################################################
#                            Copy functions
#
####################################################################################################
import Base.copy
function copy(ϕ::LatticeSite)
    LatticeSite([ϕ.A[1],ϕ.A[2],ϕ.A[3]],ϕ.θ⁺,ϕ.θ⁻,ϕ.u⁺,ϕ.u⁻)
end
function copy(c::SystConstants)
    SystConstants(c.L, c.L₃, c.g⁻², c.ν, c.κ₅, c.f, c.β)
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
function copy(ψ_list::Array{State,1})
    return [copy(ψ) for ψ in ψ_list]
end
function copy(sim::Controls)
    return Controls(sim.θmax, sim.umax, sim.Amax)
end
function copy(sim_list::Array{Controls,1})
    return [copy(sim) for sim in sim_list]
end

# -------------------------------------------------------------------------------------------------
# Comparison operator for custom types

import Base.==
function ==(ϕ₁::LatticeSite, ϕ₂::LatticeSite)
    return (ϕ₁.A[1] == ϕ₂.A[1] && ϕ₁.A[2] == ϕ₂.A[2] && ϕ₁.A[3] == ϕ₂.A[3] && ϕ₁.θ⁺ == ϕ₂.θ⁺ && 
        ϕ₁.θ⁻ == ϕ₂.θ⁻ && ϕ₁.u⁺ == ϕ₂.u⁺ && ϕ₁.u⁻ == ϕ₂.u⁻)
end
function ==(c₁::SystConstants, c₂::SystConstants)
    return (c₁.L == c₂.L && c₁.L₃ == c₂.L₃ && c₁.g⁻² == c₂.g⁻² && c₁.ν == c₂.ν
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

	consts = SystConstants(L, L, 1/0.3^2, 0.3, 1.0, 2.0/L, 0.5)
    # Construct ordered state 
    if choice == 1
        
        # Construct LxL lattice of LxL LatticeSites
        lattice = [LatticeSite([0,0,0],0,0,1.0,0.0) for y=1:L, x=1:L, z=1:L]
		nbl = latticeNeighbors(lattice,L,L)
		nnbl = latticeNextNeighbors(lattice,L,L)
		nnnbl = latticeNNNeighbors(lattice,L,L)
        ψ = State(lattice, consts, nbl, nnbl, nnnbl)
        
    # Construct random state
    elseif choice == 2
        Amax::Int64 = 2^6
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
							  rand(Uniform(0,2π)), rand(Uniform(0,2π)), 1.0, 0.0) for y=1:L, x=1:L, z=1:L]
#		lattice = [LatticeSite([0.0,0.0,0.0],
#							  rand(Uniform(0,2π)), rand(Uniform(0,2π)), 1.0, 0.0) for y=1:L, x=1:L, z=1:L]
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

# Same as above but with non-default parameter-inputs from SystConstants and options
# 1:    Correlated mean field state with specified values for fields.
# 2:    Completely random state for all fields.
# 3:    u⁺ is random, the rest are mean field.
# 4:    u⁺ and u⁻ are random, the rest are mean field.
# 5:    θ⁺ and θ⁻ are random, the rest are mean field.
function State(choice::Int64, consts::SystConstants; u⁺=1.0, u⁻=0.0, θ⁺=0.0, θ⁻=0.0, A=[0.0, 0.0, 0.0])
    L = consts.L
    L₃ = consts.L₃
    L <= 4 && throw(DomainError())
    Amax::Int64 = 10
    umax::Int64 = 3
    # Construct ordered state 
    if choice == 1
        # Mean field lattice sites for all fields.
        # Construct NxN lattice of NxN LatticeSites
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, u⁺, u⁻) for y=1:L, x=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    # Construct random state
    elseif choice == 2
		lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
                              rand(Uniform(0,2π)), rand(Uniform(0,2π)), umax*rand(), umax*rand()) for y=1:L, x=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    elseif choice == 3
        # Uniform mean field for all fields except u⁺ which is random.
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, umax*rand(), u⁻) for v=1:L, h=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    elseif choice == 4
        # Uniform mean field except for u⁺ and u⁻
        lattice = [LatticeSite([A[1], A[2], A[3]], θ⁺, θ⁻, umax*rand(), umax*rand()) for v=1:L, h=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    elseif choice == 5
        # Vary phases, all other fields are uniform mean fields.
        lattice = [LatticeSite([A[1], A[2], A[3]], rand(Uniform(0,2π)), rand(Uniform(0,2π)), u⁺, u⁻) for v=1:L, h=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    elseif choice == 6
        # All fields random except for amplitudes
        lattice = [LatticeSite([rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax)),rand(Uniform(-Amax,Amax))],
                               rand(Uniform(0,2π)), rand(Uniform(0,2π)), u⁺, u⁻) for v=1:L, h=1:L, z=1:L₃]
		nb = latticeNeighbors(lattice,L,L₃)
		nnb = latticeNextNeighbors(lattice,L,L₃)
		nnnb = latticeNNNeighbors(lattice,L,L₃)
        ψ = State(lattice, consts, nb, nnb, nnnb)
    # We only have choices 1 - 6 so far so other values for choice will give an error.
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
    @test typeof(ψ.consts.g⁻²) == typeof(ψ.consts.ν) == typeof(ψ.consts.f) == Float64
	@test isapprox(ψ.consts.L*abs(ψ.consts.f) % 1, 0.0, atol=0, rtol=1e-13)
end


# -------------------------------------------------------------------------------------------------
# Saves the state to file with specified filename. If mode is set to "a" then the state is appended
# to the file.
function save(ψ::State, filename::AbstractString; mode::AbstractString="w", visible=false)
    if mode != "w" && mode != "a"
        println("$(mode) is not a valid file-opening option.")
        return 0
    end
    if visible
        println("Saving state to file $(filename)")
    end
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    open(filename, mode) do f
        # Write line indicating start of state
        write(f, "state start\n")
        # Writing system constants on top of the file
        write(f, "$(ψ.consts.L)\n")
        write(f, "$(ψ.consts.L₃)\n")
#        write(f, "$(ψ.consts.γ)\n")
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
#        γ = parse(Float64, readline(file))
        g⁻² = parse(Float64, readline(file))
        ν = parse(Float64, readline(file))
        κ₅ = parse(Float64, readline(file))
        f = parse(Float64, readline(file))
        β = parse(Float64, readline(file))
        syst = SystConstants(L, g⁻², ν, f, β)
        
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
# Saves a list of states to a txt file without compression. This list is designed to give multiple
# different occasions of a state in the same system, i.e. all states are assumed to have the same
# constants, which is given in the beginning of the file.
function save(ψ_list::Array{State, 1}, filename::AbstractString="state_list.data"; different=false)
    L = ψ_list[1].consts.L
    L₃ = ψ_list[1].consts.L₃
    MAX_NUM_DIGITS = 10
    open(filename, "w") do f
        write(f, "state array\n")
        M = digits(length(ψ_list), 10, MAX_NUM_DIGITS)
        for i = 1:MAX_NUM_DIGITS
            write(f, "$(M[MAX_NUM_DIGITS-i+1])")
        end
        if different
            # First create file and save first state
            save(ψ_list[1], filename; mode="w")
            # Then go through the list and save remaining states to the same file by appending.
            for i = 2:length(ψ_list)
                save(ψ_list[i], filename; mode="a")
            end
        else
            write(f, "\n")
            write(f, "$(ψ_list[1].consts.L):$(ψ_list[1].consts.L₃):$(ψ_list[1].consts.g⁻²):")
            write(f, "$(ψ_list[1].consts.ν):$(ψ_list[1].consts.κ₅):$(ψ_list[1].consts.f):$(ψ_list[1].consts.β)\n")
            for ψ in ψ_list
                write(f, "state start\n")
                for z=1:L₃, h=1:L, v=1:L
                    ϕ = ψ.lattice[v,h,z]
                    write(f, "$(ϕ.A[1]):$(ϕ.A[2]):$(ϕ.A[3]):$(ϕ.θ⁺):$(ϕ.θ⁻):$(ϕ.u⁺):$(ϕ.u⁻)\n")
                end
            end
        end
    end
    return 1
end

# -------------------------------------------------------------------------------------------------
# Adds a state to a list already created by the save(::Array{State,1}, ::AbstractString; different=false) function.
function addToList(ψ::State, filename::AbstractString)
    MAX_NUM_DIGITS = 10 # Has to be the same as number used in save(::Array{State,1}, ::AbstractString)
    M₀ = 0
    # Get length of existing list and check that constants are the same as state.
    open(filename, "r") do file
        line = readline(file)
        if line != "state array"
            println("ERROR: Start of file is $(line)")
            throw(Domainerror())
        end
        M₀ = parse(Int64, readline(file))
        M₀ > 0 || throw(error("ERROR: Number of states in list is set to $(M₀)"))
        c_values = split(readline(file), ":")
        L = parse(Int64, c_values[1])
        L > 0 || throw(error("ERROR: States L is $(L)"))
        L₃ = parse(Int64, c_values[2])
        L₃ > 0 || throw(error("ERROR: States L₃ is $(L₃)"))
#        γ = parse(Float64, c_values[3])
        g⁻² = parse(Float64, c_values[3])
        ν = parse(Float64, c_values[4])
        κ₅ = parse(Float64, c_values[5])
        f = parse(Float64, c_values[6])
        β = parse(Float64, c_values[7])
        if !(ψ.consts.L == L && ψ.consts.L₃ == L₃ && ψ.consts.g⁻² == g⁻² && ψ.consts.ν == ν
                && ψ.consts.κ₅ == κ₅ && ψ.consts.f == f && ψ.consts.β == β)
            throw(error("ERROR: Input-state's constants do not match stored constants.
\tL\tL₃\tg⁻²\tν\tκ₅\tf\tβ
stored:\t$(L)\t$(L₃)\t$(g⁻²)\t$(ν)\t$(κ₅)\t$(f)\t$(β)
input:\t$(ψ.consts.L)\t$(ψ.consts.L₃)\t$(ψ.consts.g⁻²)\t$(ψ.consts.ν)\t$(ψ.consts.κ₅)\t$(ψ.consts.f)\t$(ψ.consts.β)"))
        end
    end
    # Save lattice to new state at the end
    open(filename, "a") do file
        write(file, "state start\n")
        for z=1:L₃, h=1:L, v=1:L
            ϕ = ψ.lattice[v,h,z]
            write(file, "$(ϕ.A[1]):$(ϕ.A[2]):$(ϕ.A[3]):$(ϕ.θ⁺):$(ϕ.θ⁻):$(ϕ.u⁺):$(ϕ.u⁻)\n")
        end
    end
    # Update length of list by 1.
    open(filename, "r+") do file
        line = readline(file)
        M_list = digits(M₀+1, 10, MAX_NUM_DIGITS)
        for i = 1:MAX_NUM_DIGITS
            write(file, "$(M_list[MAX_NUM_DIGITS-i+1])")
        end
        write(file, "\n")
    end
    return 1
end

# -------------------------------------------------------------------------------------------------
# Reads a list of states created by the save(::Array{State,1}, ::AbstractString) and
# addToList(::State, ::AbstractString) functions and returns a list of these states
function loadStates(filename::AbstractString)
    open(filename, "r") do file
        line = readline(file)
        if line != "state array"
            throw(DomainError())
        end
        M = parse(Int64, readline(file))
        M > 0 || throw(DomainError())
        ψ_l = Array{State, 1}(M)
        c_values = split(readline(file), ":")
        L = parse(Int64, c_values[1])
        L > 0 || throw(DomainError())
        L₃ = parse(Int64, c_values[2])
        L₃ > 0 || throw(DomainError())
#        γ = parse(Float64, c_values[3])
        g⁻² = parse(Float64, c_values[3])
        ν = parse(Float64, c_values[4])
        κ₅ = parse(Float64, c_values[5])
        f = parse(Float64, c_values[6])
        β = parse(Float64, c_values[7])
        for i = 1:M
            line = readline(file)
            if line != "state start"
                throw(DomainError())
            end
            syst = SystConstants(L,L₃,g⁻²,ν,κ₅,f,β)
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

            nbl = latticeNeighbors(lattice,L,L₃)
            nnbl = latticeNextNeighbors(lattice,L,L₃)
            nnnbl = latticeNNNeighbors(lattice,L,L₃)

            ψ_l[i] = State(lattice, syst, nbl, nnbl, nnnbl)
        end
        return ψ_l
    end
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

