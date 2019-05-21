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
    edge25_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge35_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner235_update = Array{LatticeSite, 1}(undef, 0)
    
    # First we update the part of internal points that has no influence on neighboring sub cuboids.
    for z = 2:l₃-1, y=2:l₂-1, x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    end
   
    # Planes:
    # Then we update internal points in plane 2 that are not in other planes
    x = 1
    for z = 2:l₃-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane2_update, (y,z,ϕ))
        end
    end
    # Same for Plane 3
    y = l₂
    for x = 2:l₁-1, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane3_update, (x,z,ϕ))
        end
    end
    # Same for Plane 5
    # # Plane 5 internal points, not in plane 2 or 3 (z=l₃, x≢1, y≢l₂)
    z = l₃
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane5_update, (x,y,ϕ))
        end
    end

    # Lines:
    # Intersection line between plane 2 and 3 except intersection point
    x = 1; y = l₂
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge23_update, (z,ϕ))
        end
    end
    # Intersection line between plane 2 and 5 except intersection point
    z = l₃; x = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge25_update, (y,ϕ))
        end
    end
    # Intersection line between plane 3 and 5 except intersection point
    z = l₃; y = l₂
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge35_update, (x,ϕ))
        end
    end

    # Point:
    # Intersection point between plane 2, 3 and 5
    x = 1; y = l₂; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner235_update, ϕ)
    end
       
    plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update
end
# Same as above but returns δE and number of lattice site changes
function updateInternalPointsRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane2_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane3_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane5_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge23_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge25_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge35_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner235_update = Array{LatticeSite, 1}(undef, 0)

    updates = 0
    en_diff = 0.0
    
    # First we update the part of internal points that has no influence on neighboring sub cuboids.
    for z = 2:l₃-1, y=2:l₂-1, x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
        end
    end
   
    # Planes:
    # Then we update internal points in plane 2 that are not in other planes
    x = 1
    for z = 2:l₃-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane2_update, (y,z,ϕ))
        end
    end
    # Same for Plane 3
    y = l₂
    for x = 2:l₁-1, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane3_update, (x,z,ϕ))
        end
    end
    # Same for Plane 5
    # # Plane 5 internal points, not in plane 2 or 3 (z=l₃, x≢1, y≢l₂)
    z = l₃
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane5_update, (x,y,ϕ))
        end
    end

    # Lines:
    # Intersection line between plane 2 and 3 except intersection point
    x = 1; y = l₂
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge23_update, (z,ϕ))
        end
    end
    # Intersection line between plane 2 and 5 except intersection point
    z = l₃; x = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge25_update, (y,ϕ))
        end
    end
    # Intersection line between plane 3 and 5 except intersection point
    z = l₃; y = l₂
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge35_update, (x,ϕ))
        end
    end

    # Point:
    # Intersection point between plane 2, 3 and 5
    x = 1; y = l₂; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner235_update, ϕ)
    end
       
    plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update, en_diff, updates
end



# Given that the internal points in a sub-cuboid has been updated to new values, since the border to the environment is
# plane 1, 4 and 6, it is plane 2, 5 and 3 that has potentially updated values which should be transfered to neighboring
# sub-cuboids.
function internalPointsTransfer(sc_nb::RemoteNeighbors{SubCuboid}, 
        plane2_update::Vector{Tuple{I,I,LatticeSite}}, plane3_update::Vector{Tuple{I,I,LatticeSite}}, 
        plane5_update::Vector{Tuple{I,I,LatticeSite}}, edge23_update::Vector{Tuple{I,LatticeSite}}, 
        edge25_update::Vector{Tuple{I,LatticeSite}}, edge35_update::Vector{Tuple{I,LatticeSite}},
        corner235_update::Vector{LatticeSite}) where I<:Int
    
    # SubCuboid scᵣ₋₁ is affected by updates in plane 2 through its shell[1]
    chan = sc_nb.rnᵣ₋₁
    remotecall_wait(updateShell1!, chan.where, chan, plane2_update, edge23_update, edge25_update, corner235_update)
    # SubCuboid scᵣ₊₃ is affected by updates in plane 5 through its shell[6]
    chan = sc_nb.rnᵣ₊₃
    remotecall_wait(updateShell6!, chan.where, chan, plane5_update, edge25_update, edge35_update, corner235_update)
    # SubCuboid scᵣ₊₂ is affected by updates in plane 3 through its shell[4]
    chan = sc_nb.rnᵣ₊₂
    remotecall_wait(updateShell4!, chan.where, chan, plane3_update, edge23_update, edge35_update, corner235_update)
    
    # SubCuboid scᵣ₋₁₊₂ is affected by updates in edge23 through its shell_edge14
    chan = sc_nb.rnᵣ₋₁₊₂
    remotecall_wait(updateShellEdge14!, chan.where, chan, edge23_update, corner235_update)
    # SubCuboid scᵣ₋₁₊₃ is affected by updates in edge25 through its shell_edge16
    chan = sc_nb.rnᵣ₋₁₊₃
    remotecall_wait(updateShellEdge16!, chan.where, chan, edge25_update, corner235_update)
    nothing
end


# Takes the remote reference to a sub-cuboid and execute the update protocal for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferInternalPoints!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    
    # 1. step: Update and transfer internal points
    plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update = updateInternalPoints!(sc)
    put!(chan, sc)
    internalPointsTransfer(sc_nb, plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update)
   
    nothing
end
function updateTransferInternalPointsRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    
    # 1. step: Update and transfer internal points
    (plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update,
     en_diff, updates) = updateInternalPointsRetEnUp!(sc)
    put!(chan, sc)
    internalPointsTransfer(sc_nb, plane2_update, plane3_update, plane5_update, edge23_update, edge25_update, edge35_update, corner235_update)
   
    en_diff, updates
end


# Protocol Step 2
# -------------------------------------------------------------------------------------------------

# This function updates the LatticeSites on the SubCuboid that are on planes 1 and 4, except
# for points with (x=l₁-1, y=1), (x=l₁, y=2), or points on intersection lines. This corresponds to the
# updates necessary for step 2 in the sub-cuboid update protocol. Updated points that could be
# categorized to be in plane 1 or 4 are placed in these variables and returned.
function updateIntersectionPlanes!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge24_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge13_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge45_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge15_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner245_update = Array{LatticeSite, 1}(undef, 0)
    corner135_update = Array{LatticeSite, 1}(undef, 0)

    
    # First we update the points in plane 4 (not x = 1 || z = l₃)
    y = 1
    for x = 2:l₁-2, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane4_update, (x,z,ϕ))
        end
    end
    # Update the 24 edge except corner 245
    y = 1; x = 1
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge24_update, (z,ϕ))
        end
    end
    # Update the 45 edge except corner
    z = l₃; y = 1
    for x = 2:l₁-2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge45_update, (x,ϕ))
        end
    end
    # To complete the update of plane 4 except intersection lines and z = 1,
    # we need to update corner 245
    x = 1; y = 1; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner245_update, ϕ)
    end

    # Update points in plane 1 (not y = l₂ || z = l₃)
    x = l₁
    for y = 3:l₂-1, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane1_update, (y,z,ϕ))
        end
    end
    # Update the 13 edge except corner
    x = l₁; y = l₂
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge13_update, (z,ϕ))
        end
    end
    # Update the 15 edge except corner
    x = l₁; z = l₃
    for y = 3:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge15_update, (y,ϕ))
        end
    end
    # To complete the update of plane 1 except intersection lines and z = 1,
    # we need to update corner 135
    x = l₁; y = l₂; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner135_update, ϕ)
    end

    plane1_update, plane4_update, edge13_update, edge15_update, edge24_update, edge45_update, corner135_update, corner245_update
end
# Same as above, but returns energy difference and the number of updates
function updateIntersectionPlanesRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge24_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge13_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge45_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge15_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner245_update = Array{LatticeSite, 1}(undef, 0)
    corner135_update = Array{LatticeSite, 1}(undef, 0)

    en_diff = 0.0
    updates = 0
    
    # First we update the points in plane 4 (not x = 1 || z = l₃)
    y = 1
    for x = 2:l₁-2, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane4_update, (x,z,ϕ))
        end
    end
    # Update the 24 edge except corner 245
    y = 1; x = 1
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge24_update, (z,ϕ))
        end
    end
    # Update the 45 edge except corner
    z = l₃; y = 1
    for x = 2:l₁-2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge45_update, (x,ϕ))
        end
    end
    # To complete the update of plane 4 except intersection lines and z = 1,
    # we need to update corner 245
    x = 1; y = 1; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner245_update, ϕ)
    end

    # Update points in plane 1 (not y = l₂ || z = l₃)
    x = l₁
    for y = 3:l₂-1, z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane1_update, (y,z,ϕ))
        end
    end
    # Update the 13 edge except corner
    x = l₁; y = l₂
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge13_update, (z,ϕ))
        end
    end
    # Update the 15 edge except corner
    x = l₁; z = l₃
    for y = 3:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge15_update, (y,ϕ))
        end
    end
    # To complete the update of plane 1 except intersection lines and z = 1,
    # we need to update corner 135
    x = l₁; y = l₂; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner135_update, ϕ)
    end

    plane1_update, plane4_update, edge13_update, edge15_update, edge24_update, edge45_update, corner135_update, corner245_update, en_diff, updates
end


# Suppose that points on the boundary planes that are not intersection lines and does not have x=l₁-1, y=1 or x=l₁, y=2,
# have been updated. This affects shell 1-4,6.
function intersectionPlanesTransfer(sc_nb::RemoteNeighbors{SubCuboid}, plane1_update::Vector{Tuple{I,I,LatticeSite}},
              plane4_update::Vector{Tuple{I,I,LatticeSite}}, edge13_update::Vector{Tuple{I,LatticeSite}},
              edge15_update::Vector{Tuple{I,LatticeSite}}, edge24_update::Vector{Tuple{I,LatticeSite}},
              edge45_update::Vector{Tuple{I,LatticeSite}}, corner135_update::Vector{LatticeSite}, 
              corner245_update::Vector{LatticeSite}) where I<:Int
    
    # SubCuboid scᵣ₋₂ is affected  by updates in plane 4 through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, plane4_update, edge24_update, edge45_update, corner245_update)
    # SubCuboid scᵣ₊₁ is affected by updates in plane 1 through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, plane1_update, edge15_update, edge13_update, corner135_update)
    
    # SubCuboid scᵣ₋₁ is affected by updates in the 24 edge through the y=1 part of its shell[1]
    chan = sc_nb.rnᵣ₋₁
    remotecall_wait(updateShell1!, chan.where, chan, edge24_update, corner245_update)
    # SubCuboid scᵣ₊₂ is affected by updates in the 13 edge through the x=l₁ part of its shell[4]
    chan = sc_nb.rnᵣ₊₂
    remotecall_wait(updateShell4!, chan.where, chan, edge13_update, corner135_update)
    # SubCuboid scᵣ₊₃ is affected by updates in the 45 edge and 15 edge through its shell[6]
    chan = sc_nb.rnᵣ₊₃
    remotecall_wait(updateShell6!, chan.where, chan, edge15_update, edge45_update, corner135_update, corner245_update)

    # SubCuboid scᵣ₋₂₊₃ is affected by updated in the 45 edge through its shell_edge36
    chan = sc_nb.rnᵣ₋₂₊₃
    remotecall_wait(updateShellEdge36!, chan.where, chan, edge45_update, corner245_update)
    # SubCuboid scᵣ₋₁₊₃ is affected by updates in the 245 corner through its shell_edge16
    chan = sc_nb.rnᵣ₋₁₊₃
    remotecall_wait(updateShellEdge16!, chan.where, chan, corner245_update)
    
    nothing
end


# Takes the remote reference to a sub-cuboid and execute the update protocal step 2 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferIntersectionPlanes!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
   
    # 2. step: Update and transfer intersection planes except intersection lines and points with (x=l₁-1, y=1)
    # and (x=l₁, y=2). For z = 1 the completely internal as well as 3 specific points are updated and transferred.
    (plane1_update, plane4_update, edge13_update, edge15_update,
     edge24_update, edge45_update, corner135_update, corner245_update) = updateIntersectionPlanes!(sc)
    put!(chan, sc)
    intersectionPlanesTransfer(sc_nb, plane1_update, plane4_update, edge13_update, edge15_update,
                               edge24_update, edge45_update, corner135_update, corner245_update)
    
    nothing
end
function updateTransferIntersectionPlanesRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
   
    # 2. step: Update and transfer intersection planes except intersection lines and points with (x=l₁-1, y=1)
    # and (x=l₁, y=2). For z = 1 the completely internal as well as 3 specific points are updated and transferred.
    (plane1_update, plane4_update, edge13_update, edge15_update,
     edge24_update, edge45_update, corner135_update, corner245_update, en_diff, updates) = updateIntersectionPlanesRetEnUp!(sc)
    put!(chan, sc)
    intersectionPlanesTransfer(sc_nb, plane1_update, plane4_update, edge13_update, edge15_update,
                               edge24_update, edge45_update, corner135_update, corner245_update)
    
    en_diff, updates
end


# Protocol Step 3
# Updates the intersection line between plane 4 and 1 in addition to two lines in each plane on each
# side of the intersection line.
# -------------------------------------------------------------------------------------------------

function updateIntersectionLines!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
   
    plane1_edge_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    plane4_edge_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge14_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge45_point_update = Array{LatticeSite, 1}(undef, 0)
    edge15_point_update = Array{LatticeSite, 1}(undef, 0)
    corner145_update = Array{LatticeSite, 1}(undef, 0)
    
    # Updating the 3 corner lines for 1 < z < l₃
    for z = 2:l₃-1
        x = l₁-1; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane4_edge_update, (z,ϕ))
        end
        
        x = l₁; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge14_update, (z,ϕ))
        end
        
        x = l₁; y = 2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane1_edge_update, (z,ϕ))
        end
    end

    # For z = l₃ we put updated points in special arrays
    z = l₃

    x = l₁-1; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(edge45_point_update, ϕ)
    end

    x = l₁; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner145_update, ϕ)
    end

    x = l₁; y = 2
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(edge15_point_update, ϕ)
    end

    plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update
end
# Same as above, but also return energy difference and number of updates
function updateIntersectionLinesRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
   
    plane1_edge_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    plane4_edge_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge14_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge45_point_update = Array{LatticeSite, 1}(undef, 0)
    edge15_point_update = Array{LatticeSite, 1}(undef, 0)
    corner145_update = Array{LatticeSite, 1}(undef, 0)

    en_diff = 0.0
    updates = 0
    
    # Updating the 3 corner lines for 1 < z < l₃
    for z = 2:l₃-1
        x = l₁-1; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane4_edge_update, (z,ϕ))
        end
        
        x = l₁; y = 1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge14_update, (z,ϕ))
        end
        
        x = l₁; y = 2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane1_edge_update, (z,ϕ))
        end
    end

    # For z = l₃ we put updated points in special arrays
    z = l₃

    x = l₁-1; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(edge45_point_update, ϕ)
    end

    x = l₁; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner145_update, ϕ)
    end

    x = l₁; y = 2
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(edge15_point_update, ϕ)
    end

    plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update, en_diff, updates
end


# Suppose that for all layers except z=1, we have updated the three points (x=l₁, y=1), (x=l₁-1, y=1) and (x=l₁, y=2).
# This function then sends any updated values to their respective sub-cube shells / shell-edges.
function intersectionLinesTransfer(sc_nb::RemoteNeighbors{SubCuboid}, plane1_edge_update::Vector{Tuple{I,LatticeSite}},
                  plane4_edge_update::Vector{Tuple{I,LatticeSite}}, edge14_update::Vector{Tuple{I,LatticeSite}},
                  edge15_point_update::Vector{LatticeSite},
                  edge45_point_update::Vector{LatticeSite}, corner145_update::Vector{LatticeSite}) where I<:Int
    
    # SubCuboid scᵣ₋₂ is affected by updates in plane 4 through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, plane4_edge_update, edge14_update, edge45_point_update, corner145_update)
    # SubCuboid scᵣ₊₁ is affected by updates in plane 1 through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, plane1_edge_update, edge14_update, edge15_point_update, corner145_update)
    # SubCuboid scᵣ₊₃ is affected by updates in plane 5 through its shell[6]
    chan = sc_nb.rnᵣ₊₃
    remotecall_wait(updateShell6!, chan.where, chan, edge15_point_update, edge45_point_update, corner145_update)
    
    # SubCuboid scᵣ₊₁₋₂ is affected by updates in edge 14 through its shell_edge23
    chan = sc_nb.rnᵣ₊₁₋₂
    remotecall_wait(updateShellEdge23!, chan.where, chan, edge14_update, corner145_update)
    # SubCuboid scᵣ₋₂₊₃ is affected by updates in edge 45 through its shell_edge36
    chan = sc_nb.rnᵣ₋₂₊₃
    remotecall_wait(updateShellEdge36!, chan.where, chan, edge45_point_update, corner145_update)
    
    nothing
end

# Takes the remote reference to a sub-cuboid and execute the update protocal step 3 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferIntersectionLines!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
   
    # 3. step: Update and transfer intersection lines in addition to lines (x=l₁-1, y=1, z>1), (x=l₁, y=2, z>1).
    plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update = updateIntersectionLines!(sc)
    put!(chan, sc)
    intersectionLinesTransfer(sc_nb, plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update)
    
    nothing
end
# Same as above, but returns energy difference and number of updates.
function updateTransferIntersectionLinesRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
   
    # 3. step: Update and transfer intersection lines in addition to lines (x=l₁-1, y=1, z>1), (x=l₁, y=2, z>1).
    (plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update,
     en_diff, updates) = updateIntersectionLinesRetEnUp!(sc)
    put!(chan, sc)
    intersectionLinesTransfer(sc_nb, plane1_edge_update, plane4_edge_update, edge14_update, edge15_point_update, edge45_point_update, corner145_update)
    
    en_diff, updates
end


# Protocol Step 4
# This is a pure z = 1 protocol step going through the non- plane1 and plane4 parts of plane6.
# -------------------------------------------------------------------------------------------------

function updatePlane6Internals!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge36_update = Array{Tuple{Int64,LatticeSite}, 1}(undef, 0)
    edge26_update = Array{Tuple{Int64,LatticeSite}, 1}(undef, 0)
    corner236_update = Array{LatticeSite, 1}(undef, 0)

    # Points in plane 6 that are not edges
    z = 1
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane6_update, (x,y,ϕ))
        end
    end

    # Edge 26 without corner
    x = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge26_update, (y,ϕ))
        end
    end

    # Edge 36 without corner
    y = l₂
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge36_update, (x,ϕ))
        end
    end

    # Updating corner 236
    x = 1; y = l₂
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner236_update, ϕ)
    end


    plane6_update, edge26_update, edge36_update, corner236_update
end
# Same as above, but also return energy difference and number of successful updates.
function updatePlane6InternalsRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    
    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge36_update = Array{Tuple{Int64,LatticeSite}, 1}(undef, 0)
    edge26_update = Array{Tuple{Int64,LatticeSite}, 1}(undef, 0)
    corner236_update = Array{LatticeSite, 1}(undef, 0)

    en_diff = 0.0
    updates = 0

    # Points in plane 6 that are not edges
    z = 1
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(plane6_update, (x,y,ϕ))
        end
    end

    # Edge 26 without corner
    x = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge26_update, (y,ϕ))
        end
    end

    # Edge 36 without corner
    y = l₂
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge36_update, (x,ϕ))
        end
    end

    # Updating corner 236
    x = 1; y = l₂
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner236_update, ϕ)
    end


    plane6_update, edge26_update, edge36_update, corner236_update, en_diff, updates
end


# Suppose we have updated all points in plane 6 that are not on the edges 46 or 16.
# This function sends information to the shells and shell-edges of neighboring sub-cuboids to reflect this.
function plane6InternalsTransfer(sc_nb::RemoteNeighbors{SubCuboid}, plane6_update::Vector{Tuple{I,I,LatticeSite}},
               edge26_update::Vector{Tuple{I,LatticeSite}}, edge36_update::Vector{Tuple{I,LatticeSite}},
               corner236_update::Vector{LatticeSite}) where I<:Int

    # SubCuboid scᵣ₋₃ is affected by plane 6 updates through its shell[5]
    chan = sc_nb.rnᵣ₋₃
    remotecall_wait(updateShell5!, chan.where, chan, plane6_update, edge26_update, edge36_update, corner236_update)
    # SubCuboid scᵣ₋₁ is affected by updates in edge 26 through its shell[1]
    chan = sc_nb.rnᵣ₋₁
    remotecall_wait(updateShell1Step4!, chan.where, chan, edge26_update, corner236_update)
    # SubCuboid scᵣ₊₂ is affected by updates in edge 36 through its shell[4]
    chan = sc_nb.rnᵣ₊₂
    remotecall_wait(updateShell4Step4!, chan.where, chan, edge36_update, corner236_update)

    # SubCuboid scᵣ₋₁₊₂ is affected by updates in corner 236 through its shell_edge14
    chan = sc_nb.rnᵣ₋₁₊₂
    remotecall_wait(updateShellEdge14!, chan.where, chan, corner236_update)
    # SubCuboid scᵣ₊₂₋₃ is affected by updates in edge 36 through its shell_edge45
    chan = sc_nb.rnᵣ₊₂₋₃
    remotecall_wait(updateShellEdge45!, chan.where, chan, edge36_update, corner236_update)
    
    nothing
end


# Takes the remote reference to a sub-cuboid and execute the update protocal step 4 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferPlane6Internals!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 4. step: Update and transfer points z=1, x ∈ [1:l₁-1], y ∈ [2:l₂] internal to plane6
    plane6_update, edge26_update, edge36_update, corner236_update = updatePlane6Internals!(sc)
    put!(chan, sc)
    plane6InternalsTransfer(sc_nb, plane6_update, edge26_update, edge36_update, corner236_update)
    
    nothing
end
# Same as above, but return also energy difference and number of successful updates.
function updateTransferPlane6InternalsRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 4. step: Update and transfer points z=1, x ∈ [1:l₁-1], y ∈ [2:l₂] internal to plane6
    plane6_update, edge26_update, edge36_update, corner236_update, en_diff, updates = updatePlane6InternalsRetEnUp!(sc)
    put!(chan, sc)
    plane6InternalsTransfer(sc_nb, plane6_update, edge26_update, edge36_update, corner236_update)
    
    en_diff, updates
end

# Protocol Step 5
# This is a pure z = 1 protocol step updating plane6 intersection lines except for x = l₁-1, y = 2
# and (x=l₁, y=1).
# -------------------------------------------------------------------------------------------------

function updatePlane6IntersectionLines!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    z = 1
    
    edge46_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge16_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner246_update = Array{LatticeSite, 1}(undef, 0)
    corner136_update = Array{LatticeSite, 1}(undef, 0)

    # Updating edge 46 except corner
    y = 1; z = 1
    for x = 2:l₁-2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge46_update, (x,ϕ))
        end
    end

    # Updating corner 246
    x = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner246_update, ϕ)
    end

    # Updating edge 16 except corner
    x = l₁; z = 1
    for y = 3:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge16_update, (y,ϕ))
        end
    end

    # Updating corner 136
    y = l₂
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner136_update, ϕ)
    end

    edge16_update, edge46_update, corner246_update, corner136_update
end
# Same as above, but also return the number of successful updates
function updatePlane6IntersectionLinesRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    z = 1
    
    edge46_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge16_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner246_update = Array{LatticeSite, 1}(undef, 0)
    corner136_update = Array{LatticeSite, 1}(undef, 0)

    en_diff = 0.0
    updates = 0

    # Updating edge 46 except corner
    y = 1; z = 1
    for x = 2:l₁-2
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge46_update, (x,ϕ))
        end
    end

    # Updating corner 246
    x = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner246_update, ϕ)
    end

    # Updating edge 16 except corner
    x = l₁; z = 1
    for y = 3:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        if updated
            en_diff += δE
            updates += 1
            push!(edge16_update, (y,ϕ))
        end
    end

    # Updating corner 136
    y = l₂
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner136_update, ϕ)
    end

    edge16_update, edge46_update, corner246_update, corner136_update, en_diff, updates
end


# Suppose we have updated all points in edge 46 and 16 except the 3 points (x=l₁-1, y=1), (x=l₁, y=1) and (x=l₁, y=2). 
# This function sends information to the shells and shell-edges of neighboring sub-cuboids to reflect this.
function plane6IntersectionLinesTransfer(sc_nb::RemoteNeighbors{SubCuboid}, edge16_update::Vector{Tuple{I,LatticeSite}},
               edge46_update::Vector{Tuple{I,LatticeSite}}, corner246_update::Vector{LatticeSite},
               corner136_update::Vector{LatticeSite}) where I<:Int

    # SubCuboid scᵣ₋₁ is affected by updates in corner 246 through its shell[1]
    chan = sc_nb.rnᵣ₋₁
    remotecall_wait(updateShell1!, chan.where, chan, corner246_update)
    # SubCuboid scᵣ₊₁ is affected by updates in edge 16 through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, edge16_update, corner136_update)
    # SubCuboid scᵣ₊₂ is affected by updates in corner 136 through its shell[4]
    chan = sc_nb.rnᵣ₊₂
    remotecall_wait(updateShell4!, chan.where, chan, corner136_update)
    # SubCuboid scᵣ₋₂ is affected by updates in edge 46 through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, edge46_update, corner246_update)
    # SubCuboid scᵣ₋₃ is affected by updates in plane 6 through its shell[5]
    chan = sc_nb.rnᵣ₋₃
    remotecall_wait(updateShell5!, chan.where, chan, edge16_update, edge46_update, corner246_update, corner136_update)

    # SubCuboid scᵣ₊₁₋₃ is affected by updates in edge 16 through its shell_edge25
    chan = sc_nb.rnᵣ₊₁₋₃
    remotecall_wait(updateShellEdge25!, chan.where, chan, edge16_update, corner136_update)
    # SubCuboid scᵣ₊₂₋₃ is affected by updates in corner 136 through its shell_edge45
    chan = sc_nb.rnᵣ₊₂₋₃
    remotecall_wait(updateShellEdge45!, chan.where, chan, corner136_update)

    nothing
end



# Takes the remote reference to a sub-cuboid and execute the update protocal step 5 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferPlane6IntersectionLines!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 5. step: Update and transfer intersection lines in plane 6 excepting intersection point and its two neighbors.
    edge16_update, edge46_update, corner246_update, corner136_update = updatePlane6IntersectionLines!(sc)
    put!(chan, sc)
    plane6IntersectionLinesTransfer(sc_nb, edge16_update, edge46_update, corner246_update, corner136_update)
    
    nothing
end
# Same as above, but also return the number of successful updates
function updateTransferPlane6IntersectionLinesRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 5. step: Update and transfer intersection lines in plane 6 excepting intersection point and its two neighbors.
    edge16_update, edge46_update, corner246_update, corner136_update, en_diff, updates = updatePlane6IntersectionLinesRetEnUp!(sc)
    put!(chan, sc)
    plane6IntersectionLinesTransfer(sc_nb, edge16_update, edge46_update, corner246_update, corner136_update)
    
    en_diff, updates
end


# Protocol Step 6
# This is a pure z = 1 protocol step updating the intersection point between planes 1, 4 and 6
# and its two neighbors along the intersection lines
# -------------------------------------------------------------------------------------------------

function updatePlane6IntersectionPoint!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    z = 1

    corner_update = Array{LatticeSite}(undef, 0)
    edge46_update = Array{LatticeSite}(undef, 0)
    edge16_update = Array{LatticeSite}(undef, 0)

    x = l₁-1; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(edge46_update, ϕ)
    end

    x = l₁; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner_update, ϕ)
    end
    
    x = l₁; y = 2
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(edge16_update, ϕ)
    end

    edge16_update, edge46_update, corner_update
end
# Same as above, but also return the number of successful updates
function updatePlane6IntersectionPointRetEnUp!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    z = 1

    corner_update = Array{LatticeSite}(undef, 0)
    edge46_update = Array{LatticeSite}(undef, 0)
    edge16_update = Array{LatticeSite}(undef, 0)

    en_diff = 0.0
    updates = 0

    x = l₁-1; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(edge46_update, ϕ)
    end

    x = l₁; y = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(corner_update, ϕ)
    end
    
    x = l₁; y = 2
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    updated, δE = updateLatticeSiteRetEn!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
    if updated
        en_diff += δE
        updates += 1
        push!(edge16_update, ϕ)
    end

    edge16_update, edge46_update, corner_update, en_diff, updates
end


# Suppose we have updated the three points (x=l₁-1, y=1), (x=l₁, y=1) and (x=l₁, y=2) in plane 6.
# This function sends information to the shells and shell-edges of neighboring sub-cuboids to reflect this.
function plane6IntersectionPointTransfer(sc_nb::RemoteNeighbors{SubCuboid}, edge16_point_update::Vector{LatticeSite},
                           edge46_point_update::Vector{LatticeSite}, corner146_update::Vector{LatticeSite}) where I<:Int

    # SubCuboid scᵣ₊₁ is affected by updates in edge 16 through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, edge16_point_update, corner146_update)
    # SubCuboid scᵣ₋₂ is affected by updates in edge 46 through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, edge46_point_update, corner146_update)
    # SubCuboid scᵣ₋₃ is affected by updates in plane 6 through its shell[5]
    chan = sc_nb.rnᵣ₋₃
    remotecall_wait(updateShell5!, chan.where, chan, edge16_point_update, edge46_point_update, corner146_update)

    # SubCuboid scᵣ₊₁₋₃ is affected by updates in edge 16 through its shell_edge25
    chan = sc_nb.rnᵣ₊₁₋₃
    remotecall_wait(updateShellEdge25!, chan.where, chan, edge16_point_update, corner146_update)
    # SubCuboid scᵣ₊₁₋₂ is affected by updates in corner 146 through its shell_edge23
    chan = sc_nb.rnᵣ₊₁₋₂
    remotecall_wait(updateShellEdge23!, chan.where, chan, corner146_update)

    nothing
end


# Takes the remote reference to a sub-cuboid and execute the update protocal step 6 for completely updating all
# lattice sites of it and transferring updated points to dependent shells in neighboring sub-cuboids.
function updateTransferPlane6IntersectionPoint!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 5. step: Update and transfer intersection lines in plane 6 excepting intersection point and its two neighbors.
    edge16_update, edge46_update, corner_update = updatePlane6IntersectionPoint!(sc)
    put!(chan, sc)
    plane6IntersectionPointTransfer(sc_nb, edge16_update, edge46_update, corner_update)
    
    nothing
end
function updateTransferPlane6IntersectionPointRetEnUp!(chan::RemoteChannel{Channel{SubCuboid}})
    sc = take!(chan)
    sc_nb = copy(sc.sc_nb)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂
   
    # 5. step: Update and transfer intersection lines in plane 6 excepting intersection point and its two neighbors.
    edge16_update, edge46_update, corner_update, en_diff, updates = updatePlane6IntersectionPointRetEnUp!(sc)
    put!(chan, sc)
    plane6IntersectionPointTransfer(sc_nb, edge16_update, edge46_update, corner_update)
    
    en_diff, updates
end
