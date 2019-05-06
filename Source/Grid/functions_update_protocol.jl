############################################################################################################################
#                               Implements the SubCuboid update protocol
#__________________________________________________________________________________________________________________________#
############################################################################################################################


# Protocol Step 1
# -------------------------------------------------------------------------------------------------
# Updates the internal points of a subcuboid, i.e. the points that are not in intersection planes, which corresponds to
# step 1 of sub-cuboid update protocol.
function updateInternalPoints!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane2_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane3_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane5_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge23_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    
    # First we update the part of internal points that has no influence on neighboring sub cuboids.
    for z = 2:l₃-1, y=2:l₂-1, x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    end
    
    # Then we update internal points in plane 2 (x=1)
    x = 1
    for z = 2:l₃, y = 2:l₂
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane2_update, (y,z,ϕ))
            if y == l₂
                push!(plane3_update, (x,z,ϕ))
                push!(edge23_update, (z,ϕ))
            end
            if z == l₃
                push!(plane5_update, (x,y,ϕ))
            end
        end
    end
    
    # Plane 3 internal points, not in plane 2 (y=l₂, x≢1)
    y = l₂
    for x = 2:l₁-1, z = 2:l₃
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane3_update, (x,z,ϕ))
            if z == l₃
                push!(plane5_update, (x,y,ϕ))
            end
        end
    end
    
    # Plane 5 internal points, not in plane 2 or 3 (z=l₃, x≢1, y≢l₂)
    z = l₃
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane5_update, (x,y,ϕ))
        end
    end
    plane2_update, plane3_update, plane5_update, edge23_update
end

# Takes the remote reference to a sub-cuboid and execute the update protocal for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferInternalPoints!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    # 1. step: Update and transfer internal points
    plane2_update, plane3_update, plane5_update, edge23_update = updateInternalPoints!(sc)
    put!(chan, sc)
    internalPointsTransfer(sc_nb, plane2_update, plane3_update, plane5_update, edge23_update)
   
    # 2. step: Update and transfer intersection planes except intersection lines and points with (x=l₁-1, y=1)
    # and (x=l₁, y=2)
    nothing
end


# Protocol Step 2
# -------------------------------------------------------------------------------------------------

# This function updates the LatticeSites on the SubCuboid that are on planes 1, 4 and 6, except
# for points with (x=l₁-1, y=1), (x=l₁, y=2), or points on intersection lines. This corresponds to the
# updates necessary for step 2 in the sub-cuboid update protocol. Updated points that could be
# categorized to be in plane 1, 4, 6, or edge23 are placed in these variables and returned.
function updateIntersectionPlanes!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge23_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    
    # First we update the points on plane 4 (not z = 1)
    y = 1
    for x = 1:l₁-2, z = 2:l₃
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane4_update, (x,z,ϕ))
        end
    end
    
    # Update points in plane 1 (not z = 1)
    x = l₁
    for y = 3:l₂, z = 2:l₃
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane1_update, (y,z,ϕ))
        end
    end
    
    # Finally update points in plane 6 (excepting intersection lines)
    z = 1
    for x = 1:l₁-1, y = 2:l₂
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane6_update, (x,y,ϕ))
            if x == 1 && y == l₂
                push!(edge23_update, (1,ϕ))
            end
        end
    end
    
    plane1_update, plane4_update, plane6_update, edge23_update
end

# Takes the remote reference to a sub-cuboid and execute the update protocal step 2 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferIntersectionPlanes!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 2. step: Update and transfer intersection planes except intersection lines and points with (x=l₁-1, y=1)
    # and (x=l₁, y=2)
    plane1_update, plane4_update, plane6_update, edge23_update = updateIntersectionPlanes!(sc)
    put!(chan, sc)
    intersectionPlanesTransfer(sc_nb, plane1_update, plane4_update, plane6_update, edge23_update, l₁, l₂)
    
    nothing
end


# Protocol Step 3
# -------------------------------------------------------------------------------------------------

function updateIntersectionLines!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge14_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    
    # Updating the 3 corner lines for z > 1
    for z = 2:l₃
        x = l₁-1; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane4_update, (x,z,ϕ))
        end
        
        x = l₁; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane4_update, (x,z,ϕ))
            push!(plane1_update, (y,z,ϕ))
            push!(edge14_update, (z,ϕ))
        end
        
        x = l₁; y = 2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane1_update, (y,z,ϕ))
        end
    end
    
    # Updating y = 1 line in z = 1 plane
    y = 1; z = 1
    for x = 1:l₁-2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane4_update, (x,z,ϕ))
            push!(plane6_update, (x,y,ϕ))
        end
    end
    
    # Updating x = l₁ line in z = 1 plane
    x = l₁; z = 1
    for y = 3:l₂
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
            push!(plane1_update, (y,z,ϕ))
            push!(plane6_update, (x,y,ϕ))
        end
    end
    plane1_update, plane4_update, plane6_update, edge14_update
end

# Takes the remote reference to a sub-cuboid and execute the update protocal step 3 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferIntersectionLines!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 3. step: Update and transfer intersection lines in addition to lines (x=l₁-1, y=1, z>1), (x=l₁, y=2, z>1)
    # and excepting these points in the z == 1 plane
    plane1_update, plane4_update, plane6_update, edge14_update = updateIntersectionLines!(sc)
    put!(chan, sc)
    intersectionLinesTransfer(sc_nb, plane1_update, plane4_update, plane6_update, edge14_update, l₁, l₂)
    
    nothing
end


# Protocol Step 4
# -------------------------------------------------------------------------------------------------

function updateIntersectionPoint!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge14_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    
       
    z = 1
    x = l₁-1; y = 1;
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
        push!(plane4_update, (x,z,ϕ))
        push!(plane6_update, (x,y,ϕ))
    end
    
    x = l₁; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
        push!(plane4_update, (x,z,ϕ))
        push!(plane1_update, (y,z,ϕ))
        push!(plane6_update, (x,y,ϕ))
        push!(edge14_update, (z,ϕ))
    end
    
    x = l₁; y = 2
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β) != 0.0
        push!(plane1_update, (y,z,ϕ))
        push!(plane6_update, (x,y,ϕ))
    end
    
    plane1_update, plane4_update, plane6_update, edge14_update
end

# Takes the remote reference to a sub-cuboid and execute the update protocal step 2 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferIntersectionPoint!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 4. step: Update and transfer intersection point and in addition points (x=l₁-1, y=1, z=1) and
    # (x=l₁, y=2, z=1)
    plane1_update, plane4_update, plane6_update, edge14_update = updateIntersectionPoint!(sc)
    put!(chan, sc)
    intersectionPointTransfer(sc_nb, plane1_update, plane4_update, plane6_update, edge14_update)
    
    nothing
end


