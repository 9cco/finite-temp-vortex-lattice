####################################################################################################
#                            Functions for constructing lattices of neighbor sites
#
####################################################################################################


# --------------------------------------------------------------------------------------------------
# Constructs mirror lattice of lattice, where each element consists of a Neighbor object listing all the
# nearest neighbor LatticeSites of the LatticeSite in the mirror lattice that is at the corresponding position
# of the Neighbor object.
function latticeNeighbors(lattice::Array{LatticeSite, 2}, L::Int64)
    nb = Array{Neighbors,2}(L,L)
    
    # Set neighbors of upper right corner
    ϕᵣ₊₁ = lattice[1,1]
    ϕᵣ₋₁ = lattice[1,L-1]
    ϕᵣ₊₂ = lattice[L,L]
    ϕᵣ₋₂ = lattice[2,L]
    nb[1,L] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    
    # Set neighbors of upper left corner
    ϕᵣ₊₁ = lattice[1,2]
    ϕᵣ₋₁ = lattice[1,L]
    ϕᵣ₊₂ = lattice[L,1]
    ϕᵣ₋₂ = lattice[2,1]
    nb[1,1] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    
    # Set neighbors of lower right corner
    ϕᵣ₊₁ = lattice[L,1]
    ϕᵣ₋₁ = lattice[L,L-1]
    ϕᵣ₊₂ = lattice[L-1,L]
    ϕᵣ₋₂ = lattice[1,L]
    nb[L,L] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    
    # Set neighbors of lower left corner
    ϕᵣ₊₁ = lattice[L,2]
    ϕᵣ₋₁ = lattice[L,L]
    ϕᵣ₊₂ = lattice[L-1,1]
    ϕᵣ₋₂ = lattice[1,1]
    nb[L,1] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    
    # Set neighbors of borders except for corners
    for i=2:L-1
        # Right y-boundary
        ϕᵣ₊₁ = lattice[i,1]
        ϕᵣ₋₁ = lattice[i,L-1]
        ϕᵣ₊₂ = lattice[i-1,L]
        ϕᵣ₋₂ = lattice[i+1,L]
        nb[i,L] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
        
        # Left y-boundary
        ϕᵣ₊₁ = lattice[i,2]
        ϕᵣ₋₁ = lattice[i,L]
        ϕᵣ₊₂ = lattice[i-1,1]
        ϕᵣ₋₂ = lattice[i+1,1]
        nb[i,1] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
        
        # Lower x-boundary
        ϕᵣ₊₁ = lattice[L,i+1]
        ϕᵣ₋₁ = lattice[L,i-1]
        ϕᵣ₊₂ = lattice[L-1,i]
        ϕᵣ₋₂ = lattice[1,i]
        nb[L,i] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
        
        # Upper x-boundary
        ϕᵣ₊₁ = lattice[1,i+1]
        ϕᵣ₋₁ = lattice[1,i-1]
        ϕᵣ₊₂ = lattice[L,i]
        ϕᵣ₋₂ = lattice[2,i]
        nb[1,i] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    end
    
    # Contribution from the bulk of lattice sites
    for x=2:(L-1), y=2:(L-1)
        ϕ = lattice[y,x]          # Lattice site at position r
        ϕᵣ₊₁ = lattice[y,x+1]       # Nearest neighbor at r+x
        ϕᵣ₋₁ = lattice[y,x-1]
        ϕᵣ₊₂ = lattice[y-1,x]       # Nearest neighbor at r+y
        ϕᵣ₋₂ = lattice[y+1,x]
        nb[y,x] = Neighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
    end
    nb
end

# --------------------------------------------------------------------------------------------------
# Constructs lattice of next nearest negihbors sites of the corresponding lattice LatticeSites
function latticeNextNeighbors(lattice::Array{LatticeSite,2}, L::Int64)
    nnb = Array{NextNeighbors,2}(L,L)
    
    # Contribution from upper right corner
    ϕᵣ₊₁₊₂ = lattice[L,1]
    ϕᵣ₊₁₋₂ = lattice[2,1]
    ϕᵣ₋₁₊₂ = lattice[L,L-1]
    ϕᵣ₋₁₋₂ = lattice[2,L-1]
    nnb[1,L] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    
    # Contribution from upper left corner
    ϕᵣ₊₁₊₂ = lattice[L,2]
    ϕᵣ₊₁₋₂ = lattice[2,2]
    ϕᵣ₋₁₊₂ = lattice[L,L]
    ϕᵣ₋₁₋₂ = lattice[2,L]
    nnb[1,1] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    
    # Contribution from lower right corner
    ϕᵣ₊₁₊₂ = lattice[L-1,1]
    ϕᵣ₊₁₋₂ = lattice[1,1]
    ϕᵣ₋₁₊₂ = lattice[L-1,L-1]
    ϕᵣ₋₁₋₂ = lattice[1,L-1]
    nnb[L,L] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    
    # Contribution from lower left corner
    ϕᵣ₊₁₊₂ = lattice[L-1,2]
    ϕᵣ₊₁₋₂ = lattice[1,2]
    ϕᵣ₋₁₊₂ = lattice[L-1,L]
    ϕᵣ₋₁₋₂ = lattice[1,L]
    nnb[L,1] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    
    # Set next nearest neighbors of borders except for corners
    for i = 2:(L-1)
        # Right y-boundary
        ϕᵣ₊₁₊₂ = lattice[i-1,1]
        ϕᵣ₊₁₋₂ = lattice[i+1,1]
        ϕᵣ₋₁₊₂ = lattice[i-1,L-1]
        ϕᵣ₋₁₋₂ = lattice[i+1,L-1]
        nnb[i,L] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Left y-boundary
        ϕᵣ₊₁₊₂ = lattice[i-1,2]
        ϕᵣ₊₁₋₂ = lattice[i+1,2]
        ϕᵣ₋₁₊₂ = lattice[i-1,L]
        ϕᵣ₋₁₋₂ = lattice[i+1,L]
        nnb[i,1] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Lower x-boundary
        ϕᵣ₊₁₊₂ = lattice[L-1,i+1]
        ϕᵣ₊₁₋₂ = lattice[1,i+1]
        ϕᵣ₋₁₊₂ = lattice[L-1,i-1]
        ϕᵣ₋₁₋₂ = lattice[1,i-1]
        nnb[L,i] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Upper x-boundary
        ϕᵣ₊₁₊₂ = lattice[L,i+1]
        ϕᵣ₊₁₋₂ = lattice[2,i+1]
        ϕᵣ₋₁₊₂ = lattice[L,i-1]
        ϕᵣ₋₁₋₂ = lattice[2,i-1]
        nnb[1,i] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    end
    
    # Contributions from bulk positions
    for h_pos = 2:(L-1), v_pos=2:(L-1)
        ϕᵣ₊₁₊₂ = lattice[v_pos-1,h_pos+1]
        ϕᵣ₊₁₋₂ = lattice[v_pos+1,h_pos+1]
        ϕᵣ₋₁₊₂ = lattice[v_pos-1,h_pos-1]
        ϕᵣ₋₁₋₂ = lattice[v_pos+1,h_pos-1]
        nnb[v_pos,h_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
    end
    nnb
end

# --------------------------------------------------------------------------------------------------
# Constructs lattice of next next nearest neighbor sites of the corresponding lattice LatticeSites
function latticeNNNeighbors(lattice::Array{LatticeSite, 2}, L::Int64)
    nnnb = Array{NNNeighbors,2}(L,L)
    
    # Contribution from upper right 2x2 corner
    for i = 0:1, j=0:1
        ϕᵣ₊₁₁ = lattice[1+i, 1+j]
        ϕᵣ₋₁₁ = lattice[1+i, L-3+j]
        ϕᵣ₊₂₂ = lattice[L-1+i, L-1+j]
        ϕᵣ₋₂₂ = lattice[3+i, L-1+j]
        nnnb[1+i,L-1+j] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
    end
    
    # Contribution from upper left 2x2 corner
    for i = 0:1, j=0:1
        ϕᵣ₊₁₁ = lattice[1+i, 3+j]
        ϕᵣ₋₁₁ = lattice[1+i, L-1+j]
        ϕᵣ₊₂₂ = lattice[L-1+i, 1+j]
        ϕᵣ₋₂₂ = lattice[3+i, 1+j]
        nnnb[1+i,1+j] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
    end
    
    # Contribution from lower left 2x2 corner
    for i = -1:0, j=0:1
        ϕᵣ₊₁₁ = lattice[L+i, 3+j]
        ϕᵣ₋₁₁ = lattice[L+i, L-1+j]
        ϕᵣ₊₂₂ = lattice[L-2+i, 1+j]
        ϕᵣ₋₂₂ = lattice[2+i, 1+j]
        nnnb[L+i,1+j] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
    end
    
    # Contribution from lower right 2x2 corner
    for i = -1:0, j=-1:0
        ϕᵣ₊₁₁ = lattice[L+i, 2+j]
        ϕᵣ₋₁₁ = lattice[L+i, L-2+j]
        ϕᵣ₊₂₂ = lattice[L-2+i, L+j]
        ϕᵣ₋₂₂ = lattice[2+i, L+j]
        nnnb[L+i,L+j] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
    end
    
    # Contribution from borders
    for k = 3:(L-2)
        for i = 0:1
            # Contribution from 2-thick upper x-boundary
            ϕᵣ₊₁₁ = lattice[1+i, k+2]
            ϕᵣ₋₁₁ = lattice[1+i, k-2]
            ϕᵣ₊₂₂ = lattice[L-1+i, k]
            ϕᵣ₋₂₂ = lattice[3+i, k]
            nnnb[1+i,k] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
            
            # Contribution from 2-thick left y-boundary
            ϕᵣ₊₁₁ = lattice[k, 3+i]
            ϕᵣ₋₁₁ = lattice[k, L-1+i]
            ϕᵣ₊₂₂ = lattice[k-2, 1+i]
            ϕᵣ₋₂₂ = lattice[k+2, 1+i]
            nnnb[k,1+i] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
        end
        
        for i = -1:0
            # Contribution from 2-thick lower x-boundary
            ϕᵣ₊₁₁ = lattice[L+i, k+2]
            ϕᵣ₋₁₁ = lattice[L+i, k-2]
            ϕᵣ₊₂₂ = lattice[L+i-2, k]
            ϕᵣ₋₂₂ = lattice[2+i, k]
            nnnb[L+i,k] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
            
            # Contribution from 2-thick right y-boundary
            ϕᵣ₊₁₁ = lattice[k, 2+i]
            ϕᵣ₋₁₁ = lattice[k, L-2+i]
            ϕᵣ₊₂₂ = lattice[k-2, L+i]
            ϕᵣ₋₂₂ = lattice[k+2, L+i]
            nnnb[k,L+i] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
        end
    end
    
    # Contribution from bulk
    for h_pos = 3:(L-2), v_pos = 3:(L-2)
        ϕᵣ₊₁₁ = lattice[v_pos, h_pos+2]
        ϕᵣ₋₁₁ = lattice[v_pos, h_pos-2]
        ϕᵣ₊₂₂ = lattice[v_pos-2, h_pos]
        ϕᵣ₋₂₂ = lattice[v_pos+2, h_pos]
        nnnb[v_pos,h_pos] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂)
    end
    nnnb
end

# --------------------------------------------------------------------------------------------------
# Checks that all the nearest neighbor references in ψ.nb are set correctly
# for a square lattice with periodic boundary conditions.
function checkNeighbors(ψ::State)
    
    n_c::Bool = true
    
    # Test neighbors of upper right corner
    n_c = n_c && ψ.nb[1,L].ϕᵣ₊₁ === ψ.lattice[1,1]
    n_c = n_c && ψ.nb[1,L].ϕᵣ₋₁ === ψ.lattice[1,L-1]
    n_c = n_c && ψ.nb[1,L].ϕᵣ₊₂ === ψ.lattice[L,L]
    n_c = n_c && ψ.nb[1,L].ϕᵣ₋₂ === ψ.lattice[2,L]

    
    # Test neighbors of upper left corner
    n_c = n_c && ψ.nb[1,1].ϕᵣ₊₁ === ψ.lattice[1,2]
    n_c = n_c && ψ.nb[1,1].ϕᵣ₋₁ === ψ.lattice[1,L]
    n_c = n_c && ψ.nb[1,1].ϕᵣ₊₂ === ψ.lattice[L,1]
    n_c = n_c && ψ.nb[1,1].ϕᵣ₋₂ === ψ.lattice[2,1]

    
    # Test neighbors of lower right corner
    n_c = n_c && ψ.nb[L,L].ϕᵣ₊₁ === ψ.lattice[L,1]
    n_c = n_c && ψ.nb[L,L].ϕᵣ₋₁ === ψ.lattice[L,L-1]
    n_c = n_c && ψ.nb[L,L].ϕᵣ₊₂ === ψ.lattice[L-1,L]
    n_c = n_c && ψ.nb[L,L].ϕᵣ₋₂ === ψ.lattice[1,L]

    
    # Test neighbors of lower left corner
    n_c = n_c && ψ.nb[L,1].ϕᵣ₊₁ === ψ.lattice[L,2]
    n_c = n_c && ψ.nb[L,1].ϕᵣ₋₁ === ψ.lattice[L,L]
    n_c = n_c && ψ.nb[L,1].ϕᵣ₊₂ === ψ.lattice[L-1,1]
    n_c = n_c && ψ.nb[L,1].ϕᵣ₋₂ === ψ.lattice[1,1]

    
    # Test neighbors of borders except for corners
    for i=2:L-1
        # Right y-boundary
        n_c = n_c && ψ.nb[i,L].ϕᵣ₊₁ === ψ.lattice[i,1]
        n_c = n_c && ψ.nb[i,L].ϕᵣ₋₁ === ψ.lattice[i,L-1]
        n_c = n_c && ψ.nb[i,L].ϕᵣ₊₂ === ψ.lattice[i-1,L]
        n_c = n_c && ψ.nb[i,L].ϕᵣ₋₂ === ψ.lattice[i+1,L]

        
        # Left y-boundary
        n_c = n_c && ψ.nb[i,1].ϕᵣ₊₁ === ψ.lattice[i,2]
        n_c = n_c && ψ.nb[i,1].ϕᵣ₋₁ === ψ.lattice[i,L]
        n_c = n_c && ψ.nb[i,1].ϕᵣ₊₂ === ψ.lattice[i-1,1]
        n_c = n_c && ψ.nb[i,1].ϕᵣ₋₂ === ψ.lattice[i+1,1]

        
        # Lower x-boundary
        n_c = n_c && ψ.nb[L,i].ϕᵣ₊₁ === ψ.lattice[L,i+1]
        n_c = n_c && ψ.nb[L,i].ϕᵣ₋₁ === ψ.lattice[L,i-1]
        n_c = n_c && ψ.nb[L,i].ϕᵣ₊₂ === ψ.lattice[L-1,i]
        n_c = n_c && ψ.nb[L,i].ϕᵣ₋₂ === ψ.lattice[1,i]

        
        # Upper x-boundary
        n_c = n_c && ψ.nb[1,i].ϕᵣ₊₁ === ψ.lattice[1,i+1]
        n_c = n_c && ψ.nb[1,i].ϕᵣ₋₁ === ψ.lattice[1,i-1]
        n_c = n_c && ψ.nb[1,i].ϕᵣ₊₂ === ψ.lattice[L,i]
        n_c = n_c && ψ.nb[1,i].ϕᵣ₋₂ === ψ.lattice[2,i]

    end
    
    # Contribution from the bulk of lattice sites
    for x=2:(L-1), y=2:(L-1)
        n_c = n_c && ψ.nb[y,x].ϕᵣ₊₁ === ψ.lattice[y,x+1]
        n_c = n_c && ψ.nb[y,x].ϕᵣ₋₁ === ψ.lattice[y,x-1]
        n_c = n_c && ψ.nb[y,x].ϕᵣ₊₂ === ψ.lattice[y-1,x]
        n_c = n_c && ψ.nb[y,x].ϕᵣ₋₂ === ψ.lattice[y+1,x]
    end
    return n_c
end
