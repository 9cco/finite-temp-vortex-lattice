############################################################################################################################
#                               Testing functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Checks for a specified position (by i₁, i₂ and i₃) that the remote neighbors to this position is correctly given by
# the fields in the supplied remote_neighbors variable.
function testRemoteNeighbors(cub::Cuboid, remote_neighbors::RemoteNeighbors{SubCuboid}, i₁::I, i₂::I, i₃::I) where I<:Int
    s₁ = size(cub.grid, 1); s₂ = size(cub.grid, 2); s₃ = size(cub.grid, 3)
    
    # Nearest neighbors
    if !(remote_neighbors.rnᵣ₊₁ == cub.grid[mod(i₁+1-1,s₁)+1, i₂, i₃])
        println("ERROR: rnᵣ₊₁ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₁ == cub.grid[mod(i₁-1-1,s₁)+1, i₂, i₃])
        println("ERROR: rnᵣ₋₁ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₊₂ == cub.grid[i₁, mod(i₂+1-1,s₂)+1, i₃])
        println("ERROR: rnᵣ₊₂ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₂ == cub.grid[i₁, mod(i₂-1-1,s₂)+1, i₃])
        println("ERROR: rnᵣ₋₂ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₊₃ == cub.grid[i₁, i₂, mod(i₃+1-1,s₃)+1])
        println("ERROR: rnᵣ₊₃ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₃ == cub.grid[i₁, i₂, mod(i₃-1-1,s₃)+1])
        println("ERROR: rnᵣ₋₃ wrong!")
        return false
    end
    
    # Next nearest neighbors
    if !(remote_neighbors.rnᵣ₊₁₋₂ == cub.grid[mod(i₁+1-1,s₁)+1, mod(i₂-1-1,s₂)+1, i₃])
        println("ERROR: rnᵣ₊₁₋₂ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₁₊₂ == cub.grid[mod(i₁-1-1,s₁)+1, mod(i₂+1-1,s₂)+1, i₃])
        println("ERROR: rnᵣ₋₁₊₂ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₊₁₋₃ == cub.grid[mod(i₁+1-1,s₁)+1, i₂, mod(i₃-1-1,s₃)+1])
        println("ERROR: rnᵣ₊₁₋₃ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₁₊₃ == cub.grid[mod(i₁-1-1,s₁)+1, i₂, mod(i₃+1-1,s₃)+1])
        println("ERROR: rnᵣ₋₁₊₃ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₊₂₋₃ == cub.grid[i₁, mod(i₂+1-1,s₂)+1, mod(i₃-1-1,s₃)+1])
        println("ERROR: rnᵣ₊₂₋₃ wrong!")
        return false
    end
    if !(remote_neighbors.rnᵣ₋₂₊₃ == cub.grid[i₁, mod(i₂-1-1,s₂)+1, mod(i₃+1-1,s₃)+1])
        println("ERROR: rnᵣ₋₂₊₃ wrong!")
        return false
    end
    true
end

# Takes some vector of shells and compares it with the shells that are stored in the SubCuboid.
# Reports and returns false if it finds some site that does not have equal values between the
# two.
function testShells(sc::SubCuboid, shells::Vector{Array{LatticeSite, 2}})
    for nr = 1:length(shells)
        for (i, ϕ) = enumerate(shells[nr])
            
            if !(sc.shell[nr][i] == ϕ)
                println("ERROR: in shell $(nr)")
                display(sc.shell[nr][i])
                display(ϕ)
                return false
            end
            
        end
    end
    true
end

# Does the same as above, but for edges in the SubCuboid.
function testShellEdges(sc::SubCuboid, shell_edge14::Vector{LatticeSite}, shell_edge16::Vector{LatticeSite},
        shell_edge23::Vector{LatticeSite}, shell_edge25::Vector{LatticeSite},
        shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite})
    
    for (i,ϕ) = enumerate(shell_edge14)
        if !(sc.shell_edge14[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge14")
            return false
        end
    end
    
    for (i,ϕ) = enumerate(shell_edge16)
        if !(sc.shell_edge16[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge16")
            return false
        end
    end
    
    for (i,ϕ) = enumerate(shell_edge23)
        if !(sc.shell_edge23[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge23")
            return false
        end
    end
    
    for (i,ϕ) = enumerate(shell_edge25)
        if !(sc.shell_edge25[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge25")
            return false
        end
    end
    
    for (i,ϕ) = enumerate(shell_edge36)
        if !(sc.shell_edge36[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge36")
            return false
        end
    end
    
    for (i,ϕ) = enumerate(shell_edge45)
        if !(sc.shell_edge45[i] == ϕ)
            println("ERROR: in $(i)th element of shell_edge45")
            return false
        end
    end
    true
end

# Helping function used in testShells(::Cuboid) to get remote neighbors
# to be used in testRemoteNeighbors(::Cuboid, ::RemoteNeighbors, ::I, ::I, ::I)
function getRemoteNeighbors(chan::RemoteChannel{Channel{SubCuboid}})
    sc = fetch(chan)
    sc.sc_nb
end

# Tests the shells of a cuboid by finding the shell-values of the neighboring SubCuboids and checking
# that these equal the values found in the SubCuboid being checked. Does this for each SubCuboid in the
# Cuboid.
function testShells(cub::Cuboid)
    s₁ = size(cub.grid, 1); s₂ = size(cub.grid, 2); s₃ = size(cub.grid, 3)
    
    for i₁ = 1:s₁, i₂ = 1:s₂, i₃ = 1:s₃
        
        # First we find and test remote neighbors of the subcuboid
        sc_chan = cub.grid[i₁, i₂, i₃]
        remote_neighbors = getRemoteNeighbors(sc_chan)
        if !testRemoteNeighbors(cub, remote_neighbors, i₁, i₂, i₃)
            println("ERROR: in remote neighbors in sub cuboid [$(i₁), $(i₂), $(i₃)]!")
            return false
        end
        
        # Then we use the remote neighbors to construct shells and shell edges
        shells = getShells(remote_neighbors)
        shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45 = getShellEdges(remote_neighbors)
        
        # Then we check that this matches the shells in the sub-cuboid
        sc = fetch(sc_chan)
        if !testShells(sc, shells)
            println("ERROR: in shells in sub cuboid [$(i₁), $(i₂), $(i₃)]!")
            return false
        end
        
        # Then test shell edges
        if !testShellEdges(sc, shell_edge14, shell_edge16, shell_edge23, shell_edge25, shell_edge36, shell_edge45)
            println("ERROR: in shell_edges of sub cuboid [$(i₁), $(i₂), $(i₃)]")
            return false
        end
    end
    true
end


############################################################################################################################
#                               Benchmarking functions
#__________________________________________________________________________________________________________________________#
############################################################################################################################

function meanAndErr(A::Vector{R}) where R<:Real
    mn = mean(A)
    er = std(A)/sqrt(length(A))
    (mn, er)
end

# --------------------------------------------------------------------------------------------------
# Takes in two numbers where the second is assumed to be the error. Rounds this to one significant
# digit and rounds the first value to the appropriate decimal place.
function scientificRounding(val::T, err::T; extra_digits::Int64 = 0) where T<:Real
    if err < 0.0
        println("Warning: error less than zero inserted: $(err)")
        err = abs(err)
    end
    # First we find the number of decimals needed for the first significant digit of the error
    # Round to first significant digit
    err_temp = round(err; sigdigits=1)#signif(err, 1)
    # Find the first significant digit
    st = "$(err_temp)"
    if st == "NaN"
        # Infinite error means we don't know what the value is.
        println("Warning: error was NaN")
        return signif(val, 1+extra_digits), err
    end
    if occursin(r"[e]", st)#ismatch(r"[e]", st)
        st = st[1]
    else
        st = reverse("$(err_temp)")
        st = match(r"[1-9]", st).match
    end
    sig = parse(Int64, st)
    # Now divide by this integer to get a number on the form 10^(-d) and extract d through log10
    digi_num = floor(Int64, -log10(err_temp/sig))
    # Finally use the number of digits to round the value (note that this number could be negative)
    val = round(val; digits=digi_num+extra_digits)
    return val, round(err; sigdigits=1+extra_digits)#signif(err, 1+extra_digits)
end

function benchmarkMcSweep(cub::Cuboid; samples = 2000, max_sec = 60)
    tr = @benchmark mcSweep!($cub) samples=samples seconds=max_sec evals=1
    s₁ = size(cub.grid, 1); s₂ = size(cub.grid, 2); s₃ = size(cub.grid, 3)
    L₁ = cub.syst.L₁; L₂ = cub.syst.L₂; L₃ = cub.syst.L₃
    times = tr.times*1e-9*1e6/(L₁*L₂*L₃)
    if length(times) != samples
        println("Warning: Not all samples done within $(max_sec) s")
    end
    mn, er = scientificRounding(meanAndErr(times)...)
    println("mcSweep used $(mn) ± $(er) μs pr. lattice site for $(L₁)×$(L₂)×$(L₃) lattice split ($(s₁),$(s₂),$(s₃))
giving a $(round(shellSize(cub)/(L₁*L₂*L₃); digits=2)) shell site to lattice site ratio")
end
