####################################################################################################
#                            Functions for constructing lattices of neighbor sites
#
####################################################################################################

# --------------------------------------------------------------------------------------------------
# Constructs mirror lattice of lattice, where each element consists of a Neighbor object listing all the
# nearest neighbor LatticeSites of the LatticeSite in the mirror lattice that is at the corresponding position
# of the Neighbor object.
function latticeNeighbors(lattice:Array{LatticeSite, 3}, L::Int64, L₃::Int64)
    nb = Arrray{NearestNeighbors,3}(L,L,L₃)

    for z_pos = 1:L₃, v_pos=1:L, h_pos=1:L
        ϕᵣ₊₁ = lattice[v_pos, mod(h_pos+1-1,L)+1, z_pos]
        ϕᵣ₋₁ = lattice[v_pos, mod(h_pos-1-1,L)+1, z_pos]
        ϕᵣ₊₂ = lattice[mod(v_pos-1-1,L)+1, h_pos, z_pos]
        ϕᵣ₋₂ = lattice[mod(v_pos+1-1,L)+1, h_pos, z_pos]
        ϕᵣ₊₃ = lattice[v_pos, h_pos, mod(z_pos-1-1,L₃)+1]
        ϕᵣ₋₃ = lattice[v_pos, h_pos, mod(z_pos+1-1,L₃)+1]
        nb[v_pos, h_pos, z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂, ϕᵣ₊₃, ϕᵣ₋₃)
    end
    nb
end

# --------------------------------------------------------------------------------------------------
# Old version of the above, depricated.
#function latticeNeighbors(lattice::Array{LatticeSite, 3}, L::Int64, L₃::Int64)
#    nb = Array{NearestNeighbors,3}(L,L,L₃)
#    
#    # All layers should be set in the same way.
#    for z_pos = 1:L₃
#        # Set neighbors of upper right corner
#        ϕᵣ₊₁ = lattice[1,1,z_pos]
#        ϕᵣ₋₁ = lattice[1,L-1,z_pos]
#        ϕᵣ₊₂ = lattice[L,L,z_pos]
#        ϕᵣ₋₂ = lattice[2,L,z_pos]
#        nb[1,L,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        
#        # Set neighbors of upper left corner
#        ϕᵣ₊₁ = lattice[1,2,z_pos]
#        ϕᵣ₋₁ = lattice[1,L,z_pos]
#        ϕᵣ₊₂ = lattice[L,1,z_pos]
#        ϕᵣ₋₂ = lattice[2,1,z_pos]
#        nb[1,1,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        
#        # Set neighbors of lower right corner
#        ϕᵣ₊₁ = lattice[L,1,z_pos]
#        ϕᵣ₋₁ = lattice[L,L-1,z_pos]
#        ϕᵣ₊₂ = lattice[L-1,L,z_pos]
#        ϕᵣ₋₂ = lattice[1,L,z_pos]
#        nb[L,L,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        
#        # Set neighbors of lower left corner
#        ϕᵣ₊₁ = lattice[L,2,z_pos]
#        ϕᵣ₋₁ = lattice[L,L,z_pos]
#        ϕᵣ₊₂ = lattice[L-1,1,z_pos]
#        ϕᵣ₋₂ = lattice[1,1,z_pos]
#        nb[L,1,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        
#        # Set neighbors of borders except for corners
#        for i=2:L-1
#            # Right y-boundary
#            ϕᵣ₊₁ = lattice[i,1,z_pos]
#            ϕᵣ₋₁ = lattice[i,L-1,z_pos]
#            ϕᵣ₊₂ = lattice[i-1,L,z_pos]
#            ϕᵣ₋₂ = lattice[i+1,L,z_pos]
#            nb[i,L,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#            
#            # Left y-boundary
#            ϕᵣ₊₁ = lattice[i,2,z_pos]
#            ϕᵣ₋₁ = lattice[i,L,z_pos]
#            ϕᵣ₊₂ = lattice[i-1,1,z_pos]
#            ϕᵣ₋₂ = lattice[i+1,1,z_pos]
#            nb[i,1,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#            
#            # Lower x-boundary
#            ϕᵣ₊₁ = lattice[L,i+1,z_pos]
#            ϕᵣ₋₁ = lattice[L,i-1,z_pos]
#            ϕᵣ₊₂ = lattice[L-1,i,z_pos]
#            ϕᵣ₋₂ = lattice[1,i,z_pos]
#            nb[L,i,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#            
#            # Upper x-boundary
#            ϕᵣ₊₁ = lattice[1,i+1,z_pos]
#            ϕᵣ₋₁ = lattice[1,i-1,z_pos]
#            ϕᵣ₊₂ = lattice[L,i,z_pos]
#            ϕᵣ₋₂ = lattice[2,i,z_pos]
#            nb[1,i,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        end
#        
#        # Contribution from the bulk of lattice sites
#        for x=2:(L-1), y=2:(L-1)
#            ϕ = lattice[y,x,z_pos]          # Lattice site at position r
#            ϕᵣ₊₁ = lattice[y,x+1,z_pos]       # Nearest neighbor at r+x
#            ϕᵣ₋₁ = lattice[y,x-1,z_pos]
#            ϕᵣ₊₂ = lattice[y-1,x,z_pos]       # Nearest neighbor at r+y
#            ϕᵣ₋₂ = lattice[y+1,x,z_pos]
#            nb[y,x,z_pos] = NearestNeighbors(ϕᵣ₊₁, ϕᵣ₋₁, ϕᵣ₊₂, ϕᵣ₋₂)
#        end
#    end
#    nb
#end

# --------------------------------------------------------------------------------------------------
# Constructs lattice of next nearest negihbors sites of the corresponding lattice LatticeSites
function latticeNextNeighbors(lattice::Array{LatticeSite,3}, L::Int64, L₃::Int64)
    nnb = Array{NextNeighbors,3}(L,L,L₃)
    
    for z_pos = 1:L₃
        # Contribution from upper right corner
        ϕᵣ₊₁₊₂ = lattice[L,1,z_pos]
        ϕᵣ₊₁₋₂ = lattice[2,1,z_pos]
        ϕᵣ₋₁₊₂ = lattice[L,L-1,z_pos]
        ϕᵣ₋₁₋₂ = lattice[2,L-1,z_pos]
        nnb[1,L,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Contribution from upper left corner
        ϕᵣ₊₁₊₂ = lattice[L,2,z_pos]
        ϕᵣ₊₁₋₂ = lattice[2,2,z_pos]
        ϕᵣ₋₁₊₂ = lattice[L,L,z_pos]
        ϕᵣ₋₁₋₂ = lattice[2,L,z_pos]
        nnb[1,1,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Contribution from lower right corner
        ϕᵣ₊₁₊₂ = lattice[L-1,1,z_pos]
        ϕᵣ₊₁₋₂ = lattice[1,1,z_pos]
        ϕᵣ₋₁₊₂ = lattice[L-1,L-1,z_pos]
        ϕᵣ₋₁₋₂ = lattice[1,L-1,z_pos]
        nnb[L,L,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Contribution from lower left corner
        ϕᵣ₊₁₊₂ = lattice[L-1,2,z_pos]
        ϕᵣ₊₁₋₂ = lattice[1,2,z_pos]
        ϕᵣ₋₁₊₂ = lattice[L-1,L,z_pos]
        ϕᵣ₋₁₋₂ = lattice[1,L,z_pos]
        nnb[L,1,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        
        # Set next nearest neighbors of borders except for corners
        for i = 2:(L-1)
            # Right y-boundary
            ϕᵣ₊₁₊₂ = lattice[i-1,1,z_pos]
            ϕᵣ₊₁₋₂ = lattice[i+1,1,z_pos]
            ϕᵣ₋₁₊₂ = lattice[i-1,L-1,z_pos]
            ϕᵣ₋₁₋₂ = lattice[i+1,L-1,z_pos]
            nnb[i,L,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
            
            # Left y-boundary
            ϕᵣ₊₁₊₂ = lattice[i-1,2,z_pos]
            ϕᵣ₊₁₋₂ = lattice[i+1,2,z_pos]
            ϕᵣ₋₁₊₂ = lattice[i-1,L,z_pos]
            ϕᵣ₋₁₋₂ = lattice[i+1,L,z_pos]
            nnb[i,1,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
            
            # Lower x-boundary
            ϕᵣ₊₁₊₂ = lattice[L-1,i+1,z_pos]
            ϕᵣ₊₁₋₂ = lattice[1,i+1,z_pos]
            ϕᵣ₋₁₊₂ = lattice[L-1,i-1,z_pos]
            ϕᵣ₋₁₋₂ = lattice[1,i-1,z_pos]
            nnb[L,i,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
            
            # Upper x-boundary
            ϕᵣ₊₁₊₂ = lattice[L,i+1,z_pos]
            ϕᵣ₊₁₋₂ = lattice[2,i+1,z_pos]
            ϕᵣ₋₁₊₂ = lattice[L,i-1,z_pos]
            ϕᵣ₋₁₋₂ = lattice[2,i-1,z_pos]
            nnb[1,i,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        end
        
        # Contributions from bulk positions
        for h_pos = 2:(L-1), v_pos=2:(L-1)
            ϕᵣ₊₁₊₂ = lattice[v_pos-1,h_pos+1,z_pos]
            ϕᵣ₊₁₋₂ = lattice[v_pos+1,h_pos+1,z_pos]
            ϕᵣ₋₁₊₂ = lattice[v_pos-1,h_pos-1,z_pos]
            ϕᵣ₋₁₋₂ = lattice[v_pos+1,h_pos-1,z_pos]
            nnb[v_pos,h_pos,z_pos] = NextNeighbors(ϕᵣ₊₁₊₂, ϕᵣ₊₁₋₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂)
        end
    end
    nnb
end

# --------------------------------------------------------------------------------------------------
# Constructs lattice of next next nearest neighbor sites of the corresponding lattice LatticeSites
function latticeNNNeighbors(lattice::Array{LatticeSite, 3}, L::Int64, L₃::Int64)
    nnnb = Array{NNNeighbors,3}(L,L,L₃)
    for z_pos = 1:L₃, h_pos=1:L, v_pos=1:L
        ϕᵣ₊₁₁ = lattice[v_pos, mod(h_pos+2-1,L)+1, z_pos]
        ϕᵣ₋₁₁ = lattice[v_pos, mod(h_pos-2-1,L)+1, z_pos]
        ϕᵣ₊₂₂ = lattice[mod(v_pos-2-1,L)+1, h_pos, z_pos]
        ϕᵣ₋₂₂ = lattice[mod(v_pos+2-1,L)+1, h_pos, z_pos]
        ϕᵣ₊₃₃ = lattice[v_pos, h_pos, mod(z_pos-3-1,L)+1]
        ϕᵣ₋₃₃ = lattice[v_pos, h_pos, mod(z_pos+3-1,L)+1]
        nnnb[v_pos, h_pos, z_pos] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
    end
    nnnb
end

# An earlier iteration of the above which is less time-consuming since there is no need to calculate
# the modulus, but has much more complicated code
#function latticeNNNeighbors(lattice::Array{LatticeSite, 3}, L::Int64, L₃::Int64)
#    nnnb = Array{NNNeighbors,3}(L,L,L₃)
#
#    # First we enumerate the bulk layers of the cube where next nearest neighbors in the z-direction
#    # are elementary
#    for z_pos = 3:L₃-2
#    
#        # Contribution from upper right 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 1+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-3+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, L-1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, L-1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, L-1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[1+i, L-1+j, z_pos+3]
#            nnnb[1+i,L-1+j, z_pos] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from upper left 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, 1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[1+i, 1+j, z_pos+3]
#            nnnb[1+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower left 2x2 corner
#        for i = -1:0, j=0:1
#            ϕᵣ₊₁₁ = lattice[L+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, 1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[L+i, 1+j, z_pos+3]
#            nnnb[L+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower right 2x2 corner
#        for i = -1:0, j=-1:0
#            ϕᵣ₊₁₁ = lattice[L+i, 2+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-2+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, L+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, L+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, L+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[L+i, L+j, z_pos+3]
#            nnnb[L+i,L+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from borders
#        for k = 3:(L-2)
#            for i = 0:1
#                # Contribution from 2-thick upper x-boundary
#                ϕᵣ₊₁₁ = lattice[1+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[1+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L-1+i, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[3+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[1+i, k, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[1+i, k, z_pos+3]
#                nnnb[1+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick left y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 3+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-1+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, 1+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, 1+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, 1+i, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[k, 1+i, z_pos+3]
#                nnnb[k,1+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#            
#            for i = -1:0
#                # Contribution from 2-thick lower x-boundary
#                ϕᵣ₊₁₁ = lattice[L+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[L+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L+i-2, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[2+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[L+i, k, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[L+i, k, z_pos+3]
#                nnnb[L+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick right y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 2+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-2+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, L+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, L+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, L+i, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[k, L+i, z_pos+3]
#                nnnb[k,L+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#        end
#        
#        # Contribution from bulk
#        for h_pos = 3:(L-2), v_pos = 3:(L-2)
#            ϕᵣ₊₁₁ = lattice[v_pos, h_pos+2, z_pos]
#            ϕᵣ₋₁₁ = lattice[v_pos, h_pos-2, z_pos]
#            ϕᵣ₊₂₂ = lattice[v_pos-2, h_pos, z_pos]
#            ϕᵣ₋₂₂ = lattice[v_pos+2, h_pos, z_pos]
#            ϕᵣ₊₃₃ = lattice[v_pos, h_pos, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[v_pos, h_pos, z_pos+3]
#            nnnb[v_pos,h_pos, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#    end
#
#    # Then we enumerate the top two layers
#    for z_pos = 1:2
#    
#        # Contribution from upper right 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 1+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-3+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, L-1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, L-1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, L-1+j, L₃+z_pos-2]
#            ϕᵣ₋₃₃ = lattice[1+i, L-1+j, z_pos+3]
#            nnnb[1+i,L-1+j, z_pos] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from upper left 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, 1+j, L₃+z_pos-2]
#            ϕᵣ₋₃₃ = lattice[1+i, 1+j, z_pos+3]
#            nnnb[1+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower left 2x2 corner
#        for i = -1:0, j=0:1
#            ϕᵣ₊₁₁ = lattice[L+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, 1+j, L₃+z_pos-2]
#            ϕᵣ₋₃₃ = lattice[L+i, 1+j, z_pos+3]
#            nnnb[L+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower right 2x2 corner
#        for i = -1:0, j=-1:0
#            ϕᵣ₊₁₁ = lattice[L+i, 2+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-2+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, L+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, L+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, L+j, L₃+z_pos-2]
#            ϕᵣ₋₃₃ = lattice[L+i, L+j, z_pos+3]
#            nnnb[L+i,L+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from borders
#        for k = 3:(L-2)
#            for i = 0:1
#                # Contribution from 2-thick upper x-boundary
#                ϕᵣ₊₁₁ = lattice[1+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[1+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L-1+i, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[3+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[1+i, k, L₃+z_pos-2]
#                ϕᵣ₋₃₃ = lattice[1+i, k, z_pos+3]
#                nnnb[1+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick left y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 3+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-1+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, 1+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, 1+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, 1+i, L₃+z_pos-2]
#                ϕᵣ₋₃₃ = lattice[k, 1+i, z_pos+3]
#                nnnb[k,1+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#            
#            for i = -1:0
#                # Contribution from 2-thick lower x-boundary
#                ϕᵣ₊₁₁ = lattice[L+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[L+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L+i-2, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[2+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[L+i, k, L₃+z_pos-2]
#                ϕᵣ₋₃₃ = lattice[L+i, k, z_pos+3]
#                nnnb[L+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick right y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 2+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-2+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, L+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, L+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, L+i, L₃+z_pos-2]
#                ϕᵣ₋₃₃ = lattice[k, L+i, z_pos+3]
#                nnnb[k,L+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#        end
#        
#        # Contribution from bulk
#        for h_pos = 3:(L-2), v_pos = 3:(L-2)
#            ϕᵣ₊₁₁ = lattice[v_pos, h_pos+2, z_pos]
#            ϕᵣ₋₁₁ = lattice[v_pos, h_pos-2, z_pos]
#            ϕᵣ₊₂₂ = lattice[v_pos-2, h_pos, z_pos]
#            ϕᵣ₋₂₂ = lattice[v_pos+2, h_pos, z_pos]
#            ϕᵣ₊₃₃ = lattice[v_pos, h_pos, L₃+z_pos-2]
#            ϕᵣ₋₃₃ = lattice[v_pos, h_pos, z_pos+3]
#            nnnb[v_pos,h_pos, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#    end
#
#    # Finally we enumerate the lower 2 layers
#    for z_pos = L₃-1:L₃
#    
#        # Contribution from upper right 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 1+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-3+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, L-1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, L-1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, L-1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[1+i, L-1+j, z_pos-L₃+2]
#            nnnb[1+i,L-1+j, z_pos] = NNNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from upper left 2x2 corner
#        for i = 0:1, j=0:1
#            ϕᵣ₊₁₁ = lattice[1+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[1+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-1+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[3+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[1+i, 1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[1+i, 1+j, z_pos-L₃+2]
#            nnnb[1+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower left 2x2 corner
#        for i = -1:0, j=0:1
#            ϕᵣ₊₁₁ = lattice[L+i, 3+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-1+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, 1+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, 1+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, 1+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[L+i, 1+j, z_pos-L₃+2]
#            nnnb[L+i,1+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from lower right 2x2 corner
#        for i = -1:0, j=-1:0
#            ϕᵣ₊₁₁ = lattice[L+i, 2+j, z_pos]
#            ϕᵣ₋₁₁ = lattice[L+i, L-2+j, z_pos]
#            ϕᵣ₊₂₂ = lattice[L-2+i, L+j, z_pos]
#            ϕᵣ₋₂₂ = lattice[2+i, L+j, z_pos]
#            ϕᵣ₊₃₃ = lattice[L+i, L+j, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[L+i, L+j, z_pos-L₃+2]
#            nnnb[L+i,L+j, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#        
#        # Contribution from borders
#        for k = 3:(L-2)
#            for i = 0:1
#                # Contribution from 2-thick upper x-boundary
#                ϕᵣ₊₁₁ = lattice[1+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[1+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L-1+i, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[3+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[1+i, k, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[1+i, k, z_pos-L₃+2]
#                nnnb[1+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick left y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 3+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-1+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, 1+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, 1+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, 1+i, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[k, 1+i, z_pos-L₃+2]
#                nnnb[k,1+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#            
#            for i = -1:0
#                # Contribution from 2-thick lower x-boundary
#                ϕᵣ₊₁₁ = lattice[L+i, k+2, z_pos]
#                ϕᵣ₋₁₁ = lattice[L+i, k-2, z_pos]
#                ϕᵣ₊₂₂ = lattice[L+i-2, k, z_pos]
#                ϕᵣ₋₂₂ = lattice[2+i, k, z_pos]
#                ϕᵣ₊₃₃ = lattice[L+i, k, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[L+i, k, z_pos-L₃+2]
#                nnnb[L+i,k, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#                
#                # Contribution from 2-thick right y-boundary
#                ϕᵣ₊₁₁ = lattice[k, 2+i, z_pos]
#                ϕᵣ₋₁₁ = lattice[k, L-2+i, z_pos]
#                ϕᵣ₊₂₂ = lattice[k-2, L+i, z_pos]
#                ϕᵣ₋₂₂ = lattice[k+2, L+i, z_pos]
#                ϕᵣ₊₃₃ = lattice[k, L+i, z_pos-3]
#                ϕᵣ₋₃₃ = lattice[k, L+i, z_pos-L₃+2]
#                nnnb[k,L+i, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#            end
#        end
#        
#        # Contribution from bulk
#        for h_pos = 3:(L-2), v_pos = 3:(L-2)
#            ϕᵣ₊₁₁ = lattice[v_pos, h_pos+2, z_pos]
#            ϕᵣ₋₁₁ = lattice[v_pos, h_pos-2, z_pos]
#            ϕᵣ₊₂₂ = lattice[v_pos-2, h_pos, z_pos]
#            ϕᵣ₋₂₂ = lattice[v_pos+2, h_pos, z_pos]
#            ϕᵣ₊₃₃ = lattice[v_pos, h_pos, z_pos-3]
#            ϕᵣ₋₃₃ = lattice[v_pos, h_pos, z_pos-L₃+2]
#            nnnb[v_pos,h_pos, z_pos] = NNeighbors(ϕᵣ₊₁₁, ϕᵣ₋₁₁, ϕᵣ₊₂₂, ϕᵣ₋₂₂, ϕᵣ₊₃₃, ϕᵣ₋₃₃)
#        end
#    end
#
#
#    nnnb
#end

# --------------------------------------------------------------------------------------------------
# Checks that all the nearest neighbor references in ψ.nb are set correctly
# for a square lattice with periodic boundary conditions.
function checkNeighbors(ψ::State)
    
    n_c::Bool = true
    L = ψ.consts.L
    L₃ = ψ.consts.L₃

    for z_pos = 1:L₃
    
        # Test neighbors of upper right corner
        n_c = n_c && ψ.nb[1,L,z_pos].ϕᵣ₊₁ === ψ.lattice[1,1,z_pos]
        n_c = n_c && ψ.nb[1,L,z_pos].ϕᵣ₋₁ === ψ.lattice[1,L-1,z_pos]
        n_c = n_c && ψ.nb[1,L,z_pos].ϕᵣ₊₂ === ψ.lattice[L,L,z_pos]
        n_c = n_c && ψ.nb[1,L,z_pos].ϕᵣ₋₂ === ψ.lattice[2,L,z_pos]

        
        # Test neighbors of upper left corner
        n_c = n_c && ψ.nb[1,1,z_pos].ϕᵣ₊₁ === ψ.lattice[1,2,z_pos]
        n_c = n_c && ψ.nb[1,1,z_pos].ϕᵣ₋₁ === ψ.lattice[1,L,z_pos]
        n_c = n_c && ψ.nb[1,1,z_pos].ϕᵣ₊₂ === ψ.lattice[L,1,z_pos]
        n_c = n_c && ψ.nb[1,1,z_pos].ϕᵣ₋₂ === ψ.lattice[2,1,z_pos]

        
        # Test neighbors of lower right corner
        n_c = n_c && ψ.nb[L,L,z_pos].ϕᵣ₊₁ === ψ.lattice[L,1,z_pos]
        n_c = n_c && ψ.nb[L,L,z_pos].ϕᵣ₋₁ === ψ.lattice[L,L-1,z_pos]
        n_c = n_c && ψ.nb[L,L,z_pos].ϕᵣ₊₂ === ψ.lattice[L-1,L,z_pos]
        n_c = n_c && ψ.nb[L,L,z_pos].ϕᵣ₋₂ === ψ.lattice[1,L,z_pos]

        
        # Test neighbors of lower left corner
        n_c = n_c && ψ.nb[L,1,z_pos].ϕᵣ₊₁ === ψ.lattice[L,2,z_pos]
        n_c = n_c && ψ.nb[L,1,z_pos].ϕᵣ₋₁ === ψ.lattice[L,L,z_pos]
        n_c = n_c && ψ.nb[L,1,z_pos].ϕᵣ₊₂ === ψ.lattice[L-1,1,z_pos]
        n_c = n_c && ψ.nb[L,1,z_pos].ϕᵣ₋₂ === ψ.lattice[1,1,z_pos]

        
        # Test neighbors of borders except for corners
        for i=2:L-1
            # Right y-boundary
            n_c = n_c && ψ.nb[i,L,z_pos].ϕᵣ₊₁ === ψ.lattice[i,1,z_pos]
            n_c = n_c && ψ.nb[i,L,z_pos].ϕᵣ₋₁ === ψ.lattice[i,L-1,z_pos]
            n_c = n_c && ψ.nb[i,L,z_pos].ϕᵣ₊₂ === ψ.lattice[i-1,L,z_pos]
            n_c = n_c && ψ.nb[i,L,z_pos].ϕᵣ₋₂ === ψ.lattice[i+1,L,z_pos]

            
            # Left y-boundary
            n_c = n_c && ψ.nb[i,1,z_pos].ϕᵣ₊₁ === ψ.lattice[i,2,z_pos]
            n_c = n_c && ψ.nb[i,1,z_pos].ϕᵣ₋₁ === ψ.lattice[i,L,z_pos]
            n_c = n_c && ψ.nb[i,1,z_pos].ϕᵣ₊₂ === ψ.lattice[i-1,1,z_pos]
            n_c = n_c && ψ.nb[i,1,z_pos].ϕᵣ₋₂ === ψ.lattice[i+1,1,z_pos]

            
            # Lower x-boundary
            n_c = n_c && ψ.nb[L,i,z_pos].ϕᵣ₊₁ === ψ.lattice[L,i+1,z_pos]
            n_c = n_c && ψ.nb[L,i,z_pos].ϕᵣ₋₁ === ψ.lattice[L,i-1,z_pos]
            n_c = n_c && ψ.nb[L,i,z_pos].ϕᵣ₊₂ === ψ.lattice[L-1,i,z_pos]
            n_c = n_c && ψ.nb[L,i,z_pos].ϕᵣ₋₂ === ψ.lattice[1,i,z_pos]

            
            # Upper x-boundary
            n_c = n_c && ψ.nb[1,i,z_pos].ϕᵣ₊₁ === ψ.lattice[1,i+1,z_pos]
            n_c = n_c && ψ.nb[1,i,z_pos].ϕᵣ₋₁ === ψ.lattice[1,i-1,z_pos]
            n_c = n_c && ψ.nb[1,i,z_pos].ϕᵣ₊₂ === ψ.lattice[L,i,z_pos]
            n_c = n_c && ψ.nb[1,i,z_pos].ϕᵣ₋₂ === ψ.lattice[2,i,z_pos]

        end
        
        # Contribution from the bulk of lattice sites
        for x=2:(L-1), y=2:(L-1)
            n_c = n_c && ψ.nb[y,x,z_pos].ϕᵣ₊₁ === ψ.lattice[y,x+1,z_pos]
            n_c = n_c && ψ.nb[y,x,z_pos].ϕᵣ₋₁ === ψ.lattice[y,x-1,z_pos]
            n_c = n_c && ψ.nb[y,x,z_pos].ϕᵣ₊₂ === ψ.lattice[y-1,x,z_pos]
            n_c = n_c && ψ.nb[y,x,z_pos].ϕᵣ₋₂ === ψ.lattice[y+1,x,z_pos]
        end
    end
    return n_c
end
