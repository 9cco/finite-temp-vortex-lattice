############################################################################################################################
#                               Neighbor functions // SubCuboidModule
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

# Finds the next nearest neighbors of points in the lattice that are internal to the different border planes. Saves these neighbors to
# the nnb 3d array.
function setBorderInternalPoints!(nnb::Array{NextNeighbors, 3}, lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}},
                                  l₁::I, l₂::I, l₃::I) where I<:Int

    # First we start with doing internal points in plane 1
    x = l₁
    for y = 2:l₂-1, z = 2:l₃-1
        ϕᵣ₊₁₋₂ = shell[1][y-1,z]
        ϕᵣ₊₁₊₂ = shell[1][y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell[1][y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Plane 2
    x = 1
    for y = 2:l₂-1, z = 2:l₃-1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = shell[2][y+1,z]
        ϕᵣ₋₁₋₂ = shell[2][y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[2][y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Plane 3
    y = l₂
    for x = 2:l₁-1, z = 2:l₃-1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = shell[3][x+1,z]
        ϕᵣ₋₁₊₂ = shell[3][x-1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell[3][x,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Plane 4
    y = 1
    for x = 2:l₁-1, z = 2:l₃-1
        ϕᵣ₊₁₋₂ = shell[4][x+1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = shell[4][x-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[4][x,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Plane 5
    z = l₃
    for x = 2:l₁-1, y = 2:l₂-1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[5][x-1,y]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[5][x,y-1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Plane 6
    z = 1
    for x = 2:l₁-1, y = 2:l₂-1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell[6][x+1,y]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell[6][x,y+1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    nothing
end

function setBorderEdges!(nnb::Array{NextNeighbors, 3}, lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}},
                         shell_edge13::Vector{LatticeSite}, shell_edge14::Vector{LatticeSite}, shell_edge16::Vector{LatticeSite},
                         shell_edge23::Vector{LatticeSite}, shell_edge24::Vector{LatticeSite},
                         shell_edge25::Vector{LatticeSite}, shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite},
                         l₁::I,l₂::I,l₃::I) where I<:Int

    # We start by doing the 4 edges parallel to the z-direction
    for z = 2:l₃-1
        # Edge 14
        x = l₁; y = 1
        ϕᵣ₊₁₋₂ = shell_edge14[z]
        ϕᵣ₊₁₊₂ = shell[1][y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = shell[4][x-1,z]
        ϕᵣ₊₁₋₃ = shell[1][y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[4][x,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 13
        x = l₁; y = l₂
        ϕᵣ₊₁₋₂ = shell[1][y-1,z]
        ϕᵣ₊₁₊₂ = shell_edge13[z]
        ϕᵣ₋₁₊₂ = shell[3][x-1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell[1][y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell[3][x,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 23
        x = 1; y = l₂
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = shell[3][x+1,z]
        ϕᵣ₋₁₊₂ = shell_edge23[z]
        ϕᵣ₋₁₋₂ = shell[2][y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[2][y,z+1]
        ϕᵣ₊₂₋₃ = shell[3][x,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 24
        x = 1; y = 1
        ϕᵣ₊₁₋₂ = shell[4][x+1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = shell[2][y+1,z]
        ϕᵣ₋₁₋₂ = shell_edge24[z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[2][y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[4][x,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Then we do edges parallel to the x-axis
    for x = 2:l₁-1
        # Edge 45
        y = 1; z = l₃
        ϕᵣ₊₁₋₂ = shell[4][x+1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = shell[4][x-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[5][x-1,y]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell_edge45[x]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 35
        y = l₂; z = l₃
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = shell[3][x+1,z]
        ϕᵣ₋₁₊₂ = shell[3][x-1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell[5][x-1,y]
        ϕᵣ₊₂₋₃ = shell[3][x,z-1]
        ϕᵣ₋₂₊₃ = shell[5][x,y-1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 36
        y = l₂; z = 1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = shell[3][x+1,z]
        ϕᵣ₋₁₊₂ = shell[3][x-1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell[6][x+1,y]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell_edge36[x]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 46
        y = 1; z = 1
        ϕᵣ₊₁₋₂ = shell[4][x+1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = shell[4][x-1,z]
        ϕᵣ₊₁₋₃ = shell[6][x+1,y]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell[6][x,y+1]
        ϕᵣ₋₂₊₃ = shell[4][x,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end

    # Finally we need to set the next nearest neighbors of points on edges parallel to the y-axis
    for y = 2:l₂-1
        # Edge 15
        x = l₁; z = l₃
        ϕᵣ₊₁₋₂ = shell[1][y-1,z]
        ϕᵣ₊₁₊₂ = shell[1][y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell[1][y,z-1]
        ϕᵣ₋₁₊₃ = shell[5][x-1,y]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[5][x,y-1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 25
        x = 1; z = l₃
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = shell[2][y+1,z]
        ϕᵣ₋₁₋₂ = shell[2][y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = shell_edge25[y]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = shell[5][x,y-1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 26
        x = 1; z = 1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = shell[2][y+1,z]
        ϕᵣ₋₁₋₂ = shell[2][y-1,z]
        ϕᵣ₊₁₋₃ = shell[6][x+1,y]
        ϕᵣ₋₁₊₃ = shell[2][y,z+1]
        ϕᵣ₊₂₋₃ = shell[6][x,y+1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

        # Edge 16
        x = l₁; z = 1
        ϕᵣ₊₁₋₂ = shell[1][y-1,z]
        ϕᵣ₊₁₊₂ = shell[1][y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = shell_edge16[y]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = shell[6][x,y+1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end
end

# Sets next nearest neighbors of the points on the 8 corners of the cuboid lattice.
function setBorderCorners!(nnb::Array{NextNeighbors, 3}, lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}},
                           shell_edge13::Vector{LatticeSite}, shell_edge14::Vector{LatticeSite}, shell_edge16::Vector{LatticeSite},
                           shell_edge23::Vector{LatticeSite}, shell_edge24::Vector{LatticeSite}, shell_edge25::Vector{LatticeSite},
                           shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite}, l₁::I,l₂::I,l₃::I) where I<:Int

    # Corners in z = l₃
    z = l₃

    # Corner 145
    x = l₁; y = 1
    ϕᵣ₊₁₋₂ = shell_edge14[z]
    ϕᵣ₊₁₊₂ = shell[1][y+1,z]
    ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
    ϕᵣ₋₁₋₂ = shell[4][x-1,z]
    ϕᵣ₊₁₋₃ = shell[1][y,z-1]
    ϕᵣ₋₁₊₃ = shell[5][x-1,y]
    ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
    ϕᵣ₋₂₊₃ = shell_edge45[x]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 135
    x = l₁; y = l₂
    ϕᵣ₊₁₋₂ = shell[1][y-1,z]
    ϕᵣ₊₁₊₂ = shell_edge13[z]
    ϕᵣ₋₁₊₂ = shell[3][x-1,z]
    ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
    ϕᵣ₊₁₋₃ = shell[1][y,z-1]
    ϕᵣ₋₁₊₃ = shell[5][x-1,y]
    ϕᵣ₊₂₋₃ = shell[3][x,z-1]
    ϕᵣ₋₂₊₃ = shell[5][x,y-1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 235
    x = 1; y = l₂
    ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
    ϕᵣ₊₁₊₂ = shell[3][x+1,z]
    ϕᵣ₋₁₊₂ = shell_edge23[z]
    ϕᵣ₋₁₋₂ = shell[2][y-1,z]
    ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
    ϕᵣ₋₁₊₃ = shell_edge25[y]
    ϕᵣ₊₂₋₃ = shell[3][x,z-1]
    ϕᵣ₋₂₊₃ = shell[5][x,y-1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 245
    x = 1; y = 1
    ϕᵣ₊₁₋₂ = shell[4][x+1,z]
    ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
    ϕᵣ₋₁₊₂ = shell[2][y+1,z]
    ϕᵣ₋₁₋₂ = shell_edge24[z]
    ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
    ϕᵣ₋₁₊₃ = shell_edge25[y]
    ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
    ϕᵣ₋₂₊₃ = shell_edge45[x]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corners in z = 1
    z = 1

    # Corner 146
    x = l₁; y = 1
    ϕᵣ₊₁₋₂ = shell_edge14[z]
    ϕᵣ₊₁₊₂ = shell[1][y+1,z]
    ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
    ϕᵣ₋₁₋₂ = shell[4][x-1,z]
    ϕᵣ₊₁₋₃ = shell_edge16[y]
    ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
    ϕᵣ₊₂₋₃ = shell[6][x,y+1]
    ϕᵣ₋₂₊₃ = shell[4][x,z+1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 136
    x = l₁; y = l₂
    ϕᵣ₊₁₋₂ = shell[1][y-1,z]
    ϕᵣ₊₁₊₂ = shell_edge13[z]
    ϕᵣ₋₁₊₂ = shell[3][x-1,z]
    ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
    ϕᵣ₊₁₋₃ = shell_edge16[y]
    ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
    ϕᵣ₊₂₋₃ = shell_edge36[x]
    ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 236
    x = 1; y = l₂
    ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
    ϕᵣ₊₁₊₂ = shell[3][x+1,z]
    ϕᵣ₋₁₊₂ = shell_edge23[z]
    ϕᵣ₋₁₋₂ = shell[2][y-1,z]
    ϕᵣ₊₁₋₃ = shell[6][x+1,y]
    ϕᵣ₋₁₊₃ = shell[2][y,z+1]
    ϕᵣ₊₂₋₃ = shell_edge36[x]
    ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)

    # Corner 246
    x = 1; y = 1
    ϕᵣ₊₁₋₂ = shell[4][x+1,z]
    ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
    ϕᵣ₋₁₊₂ = shell[2][y+1,z]
    ϕᵣ₋₁₋₂ = shell_edge24[z]
    ϕᵣ₊₁₋₃ = shell[6][x+1,y]
    ϕᵣ₋₁₊₃ = shell[2][y,z+1]
    ϕᵣ₊₂₋₃ = shell[6][x,y+1]
    ϕᵣ₋₂₊₃ = shell[4][x,z+1]
    nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
end

# ---------------------------------------------------------------------------------------------------
# Set the next nearest neighbors on the lattice of a sub-cuboid with shells around it.
function latticeNextNeighbors(lattice::Array{LatticeSite, 3}, shell::Vector{Array{LatticeSite, 2}},
                shell_edge13::Vector{LatticeSite}, shell_edge14::Vector{LatticeSite}, shell_edge16::Vector{LatticeSite},
                shell_edge23::Vector{LatticeSite}, shell_edge24::Vector{LatticeSite}, shell_edge25::Vector{LatticeSite},
                shell_edge36::Vector{LatticeSite}, shell_edge45::Vector{LatticeSite})

    l₁ = size(lattice, 1); l₂ = size(lattice, 2); l₃ = size(lattice, 3)
    (l₁ < 1 || l₂ < 1 || l₃ < 1) && throw(DomainError())
    nnb = Array{NextNeighbors, 3}(undef, l₁, l₂, l₃)
    
    # First we do all internal points, which are trivial and not dependent on the shell.
    for x = 2:l₁-1, y = 2:l₂-1, z = 2:l₃-1
        ϕᵣ₊₁₋₂ = lattice[x+1,y-1,z]
        ϕᵣ₊₁₊₂ = lattice[x+1,y+1,z]
        ϕᵣ₋₁₊₂ = lattice[x-1,y+1,z]
        ϕᵣ₋₁₋₂ = lattice[x-1,y-1,z]
        ϕᵣ₊₁₋₃ = lattice[x+1,y,z-1]
        ϕᵣ₋₁₊₃ = lattice[x-1,y,z+1]
        ϕᵣ₊₂₋₃ = lattice[x,y+1,z-1]
        ϕᵣ₋₂₊₃ = lattice[x,y-1,z+1]
        nnb[x,y,z] = NextNeighbors(ϕᵣ₊₁₋₂, ϕᵣ₊₁₊₂, ϕᵣ₋₁₊₂, ϕᵣ₋₁₋₂, ϕᵣ₊₁₋₃, ϕᵣ₋₁₊₃, ϕᵣ₊₂₋₃, ϕᵣ₋₂₊₃)
    end
    
    # Next we do points that are internal to the different planes
    setBorderInternalPoints!(nnb, lattice, shell, l₁,l₂,l₃)

    # Then we do points that are internal to edges of the lattice, i.e. the edges except corners.
    setBorderEdges!(nnb, lattice, shell, shell_edge13, shell_edge14, shell_edge16, shell_edge23, shell_edge24,
                             shell_edge25, shell_edge36, shell_edge45, l₁,l₂,l₃)

    # Finally we set the 8 corners of the lattice
    setBorderCorners!(nnb, lattice, shell, shell_edge13, shell_edge14, shell_edge16, shell_edge23, shell_edge24,
                             shell_edge25, shell_edge36, shell_edge45, l₁,l₂,l₃)

    nnb
end


