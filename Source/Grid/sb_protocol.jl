############################################################################################################################
#                               Alternative SubCuboid update protocol
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# First we update internal sub-cuboid points as before, but then we do all the border points in serial.

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

function updateBorderPoints!(sc::SubCuboid)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃

    plane4_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge24_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge46_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge14_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge45_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner245_update = Array{LatticeSite, 1}(undef, 0)
    corner145_update = Array{LatticeSite, 1}(undef, 0)
    corner146_update = Array{LatticeSite, 1}(undef, 0)
    corner246_update = Array{LatticeSite, 1}(undef, 0)

    plane1_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge15_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge13_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge16_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner135_update = Array{LatticeSite, 1}(undef, 0)
    corner136_update = Array{LatticeSite, 1}(undef, 0)

    plane6_update = Array{Tuple{Int64,Int64, LatticeSite}, 1}(undef, 0)
    edge26_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    edge36_update = Array{Tuple{Int64, LatticeSite}, 1}(undef, 0)
    corner236_update = Array{LatticeSite, 1}(undef, 0)

    # Plane 4:
    # INTERNAL POINTS
    # We update all points internal to plane 4
    y = 1
    for z = 2:l₃-1, x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane4_update, (x,z,ϕ))
        end
    end
    # EDGES
    # Update points along the edge 24 except the end points which are corners.
    x = 1; y = 1
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge24_update, (z,ϕ))
        end
    end
    # Update edge 46 except corners
    y = 1; z = 1
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge46_update, (x,ϕ))
        end
    end
    # Update edge 14 except corners
    x = l₁; y = 1
    for z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge14_update, (z,ϕ))
        end
    end
    # Update edge 45 except corners
    y = 1; z = l₃
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge45_update, (x,ϕ))
        end
    end
    # CORNERS
    # Update corner 245
    x = 1; y = 1; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner245_update, ϕ)
    end
    # Update corner 145
    x = l₁; y = 1; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner145_update, ϕ)
    end
    # Update corner 146
    x = l₁; y = 1; z = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner146_update, ϕ)
    end
    # Update corner 246
    x = 1; y = 1; z = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner246_update, ϕ)
    end


    # Plane 1:
    # Update plane 1 internal points
    x = l₁
    for y = 2:l₂-1; z = 2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane1_update, (y,z,ϕ))
        end
    end
    # EDGES
    # Update edge 15 except corners
    x = l₁; z = l₃
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge15_update, (y,ϕ))
        end
    end
    # Update edge 13 except corners
    x = l₁; y = l₂
    for z = 1:2:l₃-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge13_update, (z,ϕ))
        end
    end
    # Update edge 16 except corners
    x = l₁; z = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge16_update, (y,ϕ))
        end
    end
    # CORNERS
    # Update corner 135
    x = l₁; y = l₂; z = l₃
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner135_update, ϕ)
    end
    # Update corner 136
    x = l₁; y = l₂; z = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner135_update, ϕ)
    end


    # Plane 6:
    # Update plane 6 internal points
    z = 1
    for x = 2:l₁-1, y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(plane6_update, (x,y,ϕ))
        end
    end
    # EDGES
    # Update edge 26 except corners
    x = 1; z = 1
    for y = 2:l₂-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge26_update, (y,ϕ))
        end
    end
    # Update edge 36 except corners
    y = l₂; z = 1
    for x = 2:l₁-1
        ϕ = sc.lattice[x,y,z]
        nb = sc.nb[x,y,z]
        nnb = sc.nnb[x,y,z]
        if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
            push!(edge36_update, (x,ϕ))
        end
    end
    # CORNER
    # Update corner 236
    x = 1; y = l₂; z = 1
    ϕ = sc.lattice[x,y,z]
    nb = sc.nb[x,y,z]
    nnb = sc.nnb[x,y,z]
    if updateLatticeSite!(ϕ, nb, nnb, sc.syst, sc.sim, sc.β)
        push!(corner236_update, ϕ)
    end

    return (plane1_update, plane4_update, plane6_update, edge13_update, edge14_update, edge15_update, edge16_update,
            edge24_update, edge26_update, edge36_update, edge46_update, edge45_update, corner135_update, corner136_update,
            corner145_update, corner146_update, corner236_update, corner245_update, corner246_update)
end

# Given that the border points in a sub-cuboid has been updated to new values, it is plane 1, 4, and 6 that has potentially
# updated values which should be transfered to neighboring sub-cuboids.
function borderPointsTransfer(sc_nb::RemoteNeighbors{SubCuboid}, plane1_update::Vector{Tuple{I,I,LatticeSite}},
        plane4_update::Vector{Tuple{I,I,LatticeSite}}, plane6_update::Vector{Tuple{I,I,LatticeSite}},
        edge13_update::Vector{Tuple{I,LatticeSite}}, edge14_update::Vector{Tuple{I,LatticeSite}},
        edge15_update::Vector{Tuple{I,LatticeSite}}, edge16_update::Vector{Tuple{I,LatticeSite}},
        edge24_update::Vector{Tuple{I,LatticeSite}}, edge26_update::Vector{Tuple{I,LatticeSite}},
        edge36_update::Vector{Tuple{I,LatticeSite}}, edge46_update::Vector{Tuple{I,LatticeSite}},
        edge45_update::Vector{Tuple{I,LatticeSite}}, corner135_update::Vector{LatticeSite},
        corner136_update::Vector{LatticeSite}, corner145_update::Vector{LatticeSite}, corner146_update::Vector{LatticeSite},
        corner236_update::Vector{LatticeSite}, corner245_update::Vector{LatticeSite}, corner246_update::Vector{LatticeSite}) where I<:Int

    # Nearest neighbour remote sub-cuboids.
    # SubCuboid scᵣ₊₁ is affected by updates in plane 1 through its shell[2]
    chan = sc_nb.rnᵣ₊₁
    remotecall_wait(updateShell2!, chan.where, chan, plane1_update, edge13_update, edge14_update, edge15_update, edge16_update,
                    corner135_update, corner136_update, corner145_update, corner146_update)
    # SubCuboid scᵣ₋₁ is affected by update in plane 4 and 6 through its shell[1]
    chan = sc_nb.rnᵣ₋₁
    remotecall_wait(updateShell1!, chan.where, chan, edge24_update, edge26_update, corner236_update, corner245_update, corner246_update)
    # SubCuboid scᵣ₊₂ is affected by updates in plane 1 and 6 through its shell[4]
    chan = sc_nb.rnᵣ₊₂
    remotecall_wait(updateShell4!, chan.where, chan, edge13_update, edge36_update, corner135_update, corner136_update, corner236_update)
    # SubCuboid scᵣ₋₂ is affected by update in plane 4 through its shell[3]
    chan = sc_nb.rnᵣ₋₂
    remotecall_wait(updateShell3!, chan.where, chan, plane4_update, edge14_update, edge24_update, edge54_update, edge64_update, 
                    corner145_update, corner146_update, corner245_update, corner246_update)
    # SubCuboid scᵣ₊₃ is affected by updates in plane 1 and 4 through its shell[6]
    chan = sc_nb.rnᵣ₊₃
    remotecall_wait(updateShell6!, chan.where, chan, edge15_update, edge45_update, corner135_update, corner145_update, corner245_update)
    # SubCuboid scᵣ₋₃ is affected by updates in plane 6 through its shell[5]
    chan = sc_nb.rnᵣ₋₃
    remotecall_wait(updateShell5!, chan.where, chan, plane6_update, edge16_update, edge26_update, edge36_update, edge46_update,
                    corner136_update, corner146_update, corner236_update, corner246_update)
    
    # Next nearest remote sub-cuboids.
    # SubCuboid scᵣ₊₁₊₂ is affected by updates in edge13 through its shell_edge24
    chan = sc_nb.rnᵣ₊₁₊₂
    remotecall_wait(updateShellEdge24!, chan.where, chan, edge13_update, corner135_update, corner136_update)
    # SubCuboid scᵣ₊₁₋₂ is affected by updates in plane edge14 through its shell_edge23
    chan = sc_nb.rnᵣ₊₁₋₂
    remotecall_wait(updateShellEdge23!, chan.where, chan, edge14_update, corner145_update, corner146_update)
    # SubCuboid scᵣ₊₁₋₃ is affected by updates in edge16 through its shell_edge25
    chan = sc_nb.rnᵣ₊₁₋₃
    remotecall_wait(updateShellEdge25!, chan.where, chan, edge16_update, corner136_update, corner146_update)
    # SubCuboid scᵣ₋₁₋₂ is affected by updates in edge24 through its shell_edge13
    chan = sc_nb.rnᵣ₋₁₋₂
    remotecall_wait(updateShellEdge13!, chan.where, chan, edge24_update, corner245_update, corner246_update)
    # SubCuboid scᵣ₋₁₊₃ is affected by updates in corner245 through its shell_edge16
    chan = sc_nb.rnᵣ₋₁₊₃
    remotecall_wait(updateShellEdge16!, chan.where, chan, corner245_update)
    # SubCuboid scᵣ₊₂₋₃ is affected by updates in corner136 through its shell_edge45
    chan = sc_nb.rnᵣ₊₂₋₃
    remotecall_wait(updateShellEdge45!, chan.where, chan, corner136_update)
    # SubCuboid scᵣ₋₂₊₃ is affected by updates in edge45 through its shell_edge36
    chan = sc_nb.rnᵣ₋₂₊₃
    remotecall_wait(updateShellEdge36!, chan.where, chan, edge45_update, corner145_update, corner245_update)
    
    nothing
end
