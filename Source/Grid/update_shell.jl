############################################################################################################################
#                               Provides functions for updating different parts of the shell
#__________________________________________________________________________________________________________________________#
############################################################################################################################

# Basic functions for shell and edge update functionality
function updateShell!(sc::SubCuboid, shell_nr::I, plane_update::Vector{Tuple{I,I,LatticeSite}}) where I<:Int
    shell = sc.shell[shell_nr]
    for (i, j, ϕ) in plane_update
        set!(shell[i,j], ϕ)
    end
    nothing
end

function updateShell!(sc::SubCuboid, shell_nr::I, i::I, edge_update::Vector{Tuple{I,LatticeSite}}) where I<:Int
    shell = sc.shell[shell_nr]
    for (j,ϕ) in edge_update
        set!(shell[i,j], ϕ)
    end
    nothing
end

function updateShell!(sc::SubCuboid, shell_nr::I, edge_update::Vector{Tuple{I,LatticeSite}}, j::I) where I<:Int
    shell = sc.shell[shell_nr]
    for (i,ϕ) in edge_update
        set!(shell[i,j], ϕ)
    end
    nothing
end

function updateShell!(sc::SubCuboid, shell_nr::I, corner_update::LatticeSite, i::I, j::I) where I<:Int
    shell = sc.shell[shell_nr]
    set!(shell[i,j], corner_update)
    nothing
end

function updateShell1!(chan::RemoteChannel{Channel{SubCuboid}}, corner246_update::Vector{LatticeSite})
    if length(corner246_update) > 0
        sc = take!(chan)
        l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
        shell_nr = 1

        updateShell!(sc, shell_nr, corner246_update[1], 1,1)

        put!(chan, sc)
    end
    nothing
end

function updateShell1Step5!(chan::RemoteChannel{Channel{SubCuboid}}, edge26_update::Vector{Tuple{I,LatticeSite}},
                            corner236_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 1

    updateShell!(sc, shell_nr, edge26_update, 1)
    if length(corner236_update) > 0
        updateShell!(sc, shell_nr, corner236_update[1], l₂, 1)
    end

    put!(chan, sc)
    nothing
end

function updateShell1!(chan::RemoteChannel{Channel{SubCuboid}}, edge24_update::Vector{Tuple{I,LatticeSite}},
                       corner245_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 1

    updateShell!(sc, shell_nr, 1, edge24_update)
    if length(corner245_update) > 0
        updateShell!(sc, shell_nr, corner245_update[1], 1, l₃)
    end

    put!(chan, sc)
    nothing
end

# Shell update functions tied to the update protocol.
function updateShell1!(chan::RemoteChannel{Channel{SubCuboid}}, plane2_update::Vector{Tuple{I,I,LatticeSite}},
                       edge23_update::Vector{Tuple{I,LatticeSite}}, edge25_update::Vector{Tuple{I,LatticeSite}},
                       corner235_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    updateShell!(sc, 1, plane2_update)
    updateShell!(sc, 1, sc.consts.l₂, edge23_update)
    updateShell!(sc, 1, edge25_update, l₃)
    if length(corner235_update) > 0
        updateShell!(sc, 1, corner235_update[1], l₂,l₃)
    end
    put!(chan, sc)
    nothing
end

function updateShell2!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_point_update::Vector{LatticeSite},
                       corner146_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 2

    if length(edge16_point_update) > 0
        updateShell!(sc, shell_nr, edge16_point_update[1], 2, 1)
    end
    if length(corner146_update) > 0
        updateShell!(sc, shell_nr, corner146_update[1], 1,1)
    end

    put!(chan, sc)
    nothing
end

function updateShell2!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_update::Vector{Tuple{I,LatticeSite}},
                       corner136_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 2

    updateShell!(sc, shell_nr, edge16_update, 1)
    if length(corner136_update) > 0
        updateShell!(sc, shell_nr, corner136_update[1], l₂, 1)
    end

    put!(chan, sc)
    nothing
end

function updateShell2!(chan::RemoteChannel{Channel{SubCuboid}}, plane1_edge_update::Vector{Tuple{I,LatticeSite}},
                       edge14_update::Vector{Tuple{I,LatticeSite}}, edge15_point_update::Vector{LatticeSite},
                       corner145_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 2

    updateShell!(sc, shell_nr, 2, plane1_edge_update)
    updateShell!(sc, shell_nr, 1, edge14_update)
    if length(edge15_point_update) > 0
        updateShell!(sc, shell_nr, edge15_point_update[1], 2, l₃)
    end
    if length(corner145_update) > 0
        updateShell!(sc, shell_nr, corner145_update[1], 1, l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShell2!(chan::RemoteChannel{Channel{SubCuboid}}, plane1_update::Vector{Tuple{I,I,LatticeSite}},
                       edge15_update::Vector{Tuple{I,LatticeSite}}, edge13_update::Vector{Tuple{I,LatticeSite}},
                       corner135_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 2

    updateShell!(sc, shell_nr, plane1_update)
    updateShell!(sc, shell_nr, edge15_update, l₃)
    updateShell!(sc, shell_nr, l₂, edge13_update)
    if length(corner135_update) > 0
        updateShell!(sc, shell_nr, corner135_update[1], l₂, l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShell3!(chan::RemoteChannel{Channel{SubCuboid}}, edge46_point_update::Vector{LatticeSite},
                       corner146_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 3

    if length(edge46_point_update) > 0
        updateShell!(sc, shell_nr, edge46_point_update[1], l₁-1, 1)
    end
    if length(corner146_update) > 0
        updateShell!(sc, shell_nr, corner146_update[1], l₁,1)
    end

    put!(chan, sc)
    nothing
end

function updateShell3!(chan::RemoteChannel{Channel{SubCuboid}}, edge46_update::Vector{Tuple{I,LatticeSite}},
                       corner246_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 3

    updateShell!(sc, shell_nr, edge46_update, 1)
    if length(corner246_update) > 0
        updateShell!(sc, shell_nr, corner246_update[1], 1, 1)
    end

    put!(chan, sc)
    nothing
end

function updateShell3!(chan::RemoteChannel{Channel{SubCuboid}}, plane4_edge_update::Vector{Tuple{I,LatticeSite}},
                       edge14_update::Vector{Tuple{I,LatticeSite}}, edge45_point_update::Vector{LatticeSite},
                       corner145_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 3

    updateShell!(sc, shell_nr, l₁-1, plane4_edge_update)
    updateShell!(sc, shell_nr, l₁, edge14_update)
    if length(edge45_point_update) > 0
        updateShell!(sc, shell_nr, edge45_point_update[1], l₁-1, l₃)
    end
    if length(corner145_update) > 0
        updateShell!(sc, shell_nr, corner145_update[1], l₁, l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShell3!(chan::RemoteChannel{Channel{SubCuboid}}, plane4_update::Vector{Tuple{I,I,LatticeSite}},
                       edge24_update::Vector{Tuple{I,LatticeSite}}, edge45_update::Vector{Tuple{I,LatticeSite}},
                       corner245_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 3

    updateShell!(sc, shell_nr, plane4_update)
    updateShell!(sc, shell_nr, 1, edge24_update)
    updateShell!(sc, shell_nr, edge45_update, l₃)
    if length(corner245_update) > 0
        updateShell!(sc, shell_nr, corner245_update[1], 1, l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShell4!(chan::RemoteChannel{Channel{SubCuboid}}, corner136_update::Vector{LatticeSite})
    if length(corner136_update) > 0
        sc = take!(chan)
        l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
        shell_nr = 4

        updateShell!(sc, shell_nr, corner136_update[1], l₁, 1)

        put!(chan, sc)
    end
    nothing
end

function updateShell4Step5!(chan::RemoteChannel{Channel{SubCuboid}}, edge36_update::Vector{Tuple{I,LatticeSite}},
                            corner236_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 4

    updateShell!(sc, shell_nr, edge36_update, 1)
    if length(corner236_update) > 0
        updateShell!(sc, shell_nr, corner236_update[1], 1,1)
    end

    put!(chan, sc)
    nothing
end

function updateShell4!(chan::RemoteChannel{Channel{SubCuboid}}, edge13_update::Vector{Tuple{I,LatticeSite}},
                       corner135_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 4

    updateShell!(sc, shell_nr, l₁, edge13_update)
    if length(corner135_update) > 0
        updateShell!(sc, shell_nr, corner135_update[1], l₁, l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShell4!(chan::RemoteChannel{Channel{SubCuboid}}, plane3_update::Vector{Tuple{I,I,LatticeSite}},
                       edge23_update::Vector{Tuple{I,LatticeSite}}, edge35_update::Vector{Tuple{I,LatticeSite}},
                       corner235_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₃ = sc.consts.l₃
    shell_nr = 4
    updateShell!(sc, shell_nr, plane3_update)
    updateShell!(sc, shell_nr, 1, edge23_update)
    updateShell!(sc, shell_nr, edge35_update, l₃)
    if length(corner235_update) > 0
        updateShell!(sc, shell_nr, corner235_update[1], 1, l₃)
    end
    put!(chan, sc)
    nothing
end

function updateShell5!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_point_update::Vector{LatticeSite},
                       edge46_point_update::Vector{LatticeSite}, corner146_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 5

    if length(edge16_point_update) > 0
        updateShell!(sc, shell_nr, edge16_point_update[1], l₁,2)
    end
    if length(edge46_point_update) > 0
        updateShell!(sc, shell_nr, edge46_point_update[1], l₁-1,1)
    end
    if length(corner146_update) > 0
        updateShell!(sc, shell_nr, corner146_update[1], l₁,1)
    end

    put!(chan, sc)
    nothing
end

function updateShell5Step6!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_update::Vector{Tuple{I,LatticeSite}},
                       corner136_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 5

    updateShell!(sc, shell_nr, l₁, edge16_update)
    if length(corner136_update) > 0
        updateShell!(sc, shell_nr, corner136_update[1], l₁,l₂)
    end

    put!(chan, sc)
    nothing
end

function updateShell5Step7!(chan::RemoteChannel{Channel{SubCuboid}}, edge46_update::Vector{Tuple{I,LatticeSite}},
                            corner246_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 5

    updateShell!(sc, shell_nr, edge46_update, 1)
    if length(corner246_update) > 0
        updateShell!(sc, shell_nr, corner246_update[1], 1,1)
    end

    put!(chan, sc)
    nothing
end

function updateShell5!(chan::RemoteChannel{Channel{SubCuboid}}, plane6_update::Vector{Tuple{I,I,LatticeSite}}, 
                       edge26_update::Vector{Tuple{I,LatticeSite}}, edge36_update::Vector{Tuple{I,LatticeSite}},
                       corner236_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 5

    updateShell!(sc, shell_nr, plane6_update)
    updateShell!(sc, shell_nr, 1, edge26_update)
    updateShell!(sc, shell_nr, edge36_update, l₂)
    if length(corner236_update) > 0
        updateShell!(sc, shell_nr, corner236_update[1], 1, l₂)
    end

    put!(chan, sc)
    nothing
end

function updateShell6!(chan::RemoteChannel{Channel{SubCuboid}}, edge15_point_update::Vector{LatticeSite},
                       edge45_point_update::Vector{LatticeSite}, corner145_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 6

    if length(edge15_point_update) > 0
        updateShell!(sc, shell_nr, edge15_point_update[1], l₁, 2)
    end
    if length(edge45_point_update) > 0
        updateShell!(sc, shell_nr, edge45_point_update[1], l₁-1, 1)
    end
    if length(corner145_update) > 0
        updateShell!(sc, shell_nr, corner145_update[1], l₁, 1)
    end

    put!(chan, sc)
    nothing
end

function updateShell6Step2!(chan::RemoteChannel{Channel{SubCuboid}}, edge15_update::Vector{Tuple{I,LatticeSite}},
                       corner135_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 6

    updateShell!(sc, shell_nr, l₁, edge15_update)
    if length(corner135_update) > 0
        updateShell!(sc, shell_nr, corner135_update[1], l₁, l₂)
    end

    put!(chan, sc)
    nothing
end

function updateShell6Step3!(chan::RemoteChannel{Channel{SubCuboid}}, edge45_update::Vector{Tuple{I,LatticeSite}},
                       corner245_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_nr = 6

    updateShell!(sc, shell_nr, edge45_update, 1)
    if length(corner245_update) > 0
        updateShell!(sc, shell_nr, corner245_update[1], 1, 1)
    end

    put!(chan, sc)
    nothing
end


function updateShell6!(chan::RemoteChannel{Channel{SubCuboid}}, plane5_update::Vector{Tuple{I,I,LatticeSite}},
                       edge25_update::Vector{Tuple{I,LatticeSite}}, edge35_update::Vector{Tuple{I,LatticeSite}},
                       corner235_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₂ = sc.consts.l₂
    shell_nr = 6
    updateShell!(sc, shell_nr, plane5_update)
    updateShell!(sc, shell_nr, 1, edge25_update)
    updateShell!(sc, shell_nr, edge35_update, l₂)
    if length(corner235_update) > 0
        updateShell!(sc, shell_nr, corner235_update[1], 1, l₂)
    end
    put!(chan, sc)
    nothing
end

function updateShellEdge!(shell_edge::Vector{LatticeSite}, edge_update::Vector{Tuple{I,LatticeSite}}) where I<:Int
    for (i,ϕ) in edge_update
        set!(shell_edge[i], ϕ)
    end
    nothing
end
function updateShellEdge!(shell_edge::Vector{LatticeSite}, ϕ::LatticeSite, i::I) where I<:Int
    set!(shell_edge[i], ϕ)
end

function updateShellEdge13!(chan::RemoteChannel{Channel{SubCuboid}}, corner246_update::Vector{LatticeSite})

    if length(corner246_update) > 0
        sc = take!(chan)
        shell_edge = sc.shell_edge13

        updateShellEdge!(shell_edge, corner246_update[1], 1)

        put!(chan, sc)
    end
    nothing
end

function updateShellEdge13!(chan::RemoteChannel{Channel{SubCuboid}}, edge24_update::Vector{Tuple{I,LatticeSite}},
                            corner245_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge13

    updateShellEdge!(shell_edge, edge24_update)
    if length(corner245_update) > 0
        updateShellEdge!(shell_edge, corner245_update[1], l₃)
    end
    put!(chan, sc)
    nothing
end

function updateShellEdge14!(chan::RemoteChannel{Channel{SubCuboid}}, corner236_update::Vector{LatticeSite})

    if length(corner236_update) > 0
        sc = take!(chan)
        l₃ = sc.consts.l₃
        shell_edge = sc.shell_edge14

        updateShellEdge!(shell_edge, corner236_update[1], 1)

        put!(chan, sc)
    end
    nothing
end

function updateShellEdge14!(chan::RemoteChannel{Channel{SubCuboid}}, edge23_update::Vector{Tuple{I,LatticeSite}},
                            corner235_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge14

    updateShellEdge!(shell_edge, edge23_update)
    if length(corner235_update) > 0
        updateShellEdge!(shell_edge, corner235_update[1], l₃)
    end
    put!(chan, sc)
    nothing
end

function updateShellEdge16!(chan::RemoteChannel{Channel{SubCuboid}}, corner245_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge16

    if length(corner245_update) > 0
        updateShellEdge!(shell_edge, corner245_update[1], 1)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge16!(chan::RemoteChannel{Channel{SubCuboid}}, edge25_update::Vector{Tuple{I,LatticeSite}},
                            corner235_update::Vector{LatticeSite}) where I<:Int

    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge16

    updateShellEdge!(shell_edge, edge25_update)
    if length(corner235_update) > 0
        updateShellEdge!(shell_edge, corner235_update[1], l₂)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge23!(chan::RemoteChannel{Channel{SubCuboid}}, corner146_update::Vector{LatticeSite})
    if length(corner146_update) > 0
        sc = take!(chan)
        l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
        shell_edge = sc.shell_edge23

        updateShellEdge!(shell_edge, corner146_update[1], 1)

        put!(chan, sc)
    end
    nothing
end

function updateShellEdge23!(chan::RemoteChannel{Channel{SubCuboid}}, edge14_update::Vector{Tuple{I,LatticeSite}},
                            corner145_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge23

    updateShellEdge!(shell_edge, edge14_update)
    if length(corner145_update) > 0
        updateShellEdge!(shell_edge, corner145_update[1], l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge24!(chan::RemoteChannel{Channel{SubCuboid}}, corner136_update::Vector{LatticeSite})
    if length(corner136_update) > 0
        sc = take!(chan)
        l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
        shell_edge = sc.shell_edge24

        updateShellEdge!(shell_edge, corner136_update[1], 1)

        put!(chan, sc)
    end
    nothing
end

function updateShellEdge24!(chan::RemoteChannel{Channel{SubCuboid}}, edge13_update::Vector{Tuple{I,LatticeSite}},
                            corner135_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge24

    updateShellEdge!(shell_edge, edge13_update)
    if length(corner135_update) > 0
        updateShellEdge!(shell_edge, corner135_update[1], l₃)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge25!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_point_update::Vector{LatticeSite},
                            corner146_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge25

    if length(edge16_point_update) > 0
        updateShellEdge!(shell_edge, edge16_point_update[1], 2)
    end
    if length(corner146_update) > 0
        updateShellEdge!(shell_edge, corner146_update[1], 1)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge25!(chan::RemoteChannel{Channel{SubCuboid}}, edge16_update::Vector{Tuple{I,LatticeSite}},
                            corner136_update::Vector{LatticeSite}) where I<:Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge25

    updateShellEdge!(shell_edge, edge16_update)
    if length(corner136_update) > 0
        updateShellEdge!(shell_edge, corner136_update[1], l₂)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge36!(chan::RemoteChannel{Channel{SubCuboid}}, edge45_point_update::Vector{LatticeSite},
                            corner145_update::Vector{LatticeSite})
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge36

    if length(edge45_point_update) > 0
        updateShellEdge!(shell_edge, edge45_point_update[1], l₁-1)
    end
    if length(corner145_update) > 0
        updateShellEdge!(shell_edge, corner145_update[1], l₁)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge36!(chan::RemoteChannel{Channel{SubCuboid}}, edge45_update::Vector{Tuple{I,LatticeSite}},
                            corner245_update::Vector{LatticeSite}) where I <: Int
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge36

    updateShellEdge!(shell_edge, edge45_update)
    if length(corner245_update) > 0
        updateShellEdge!(shell_edge, corner245_update[1], 1)
    end

    put!(chan, sc)
    nothing
end

function updateShellEdge45!(chan::RemoteChannel{Channel{SubCuboid}}, corner136_update::Vector{LatticeSite})
    if length(corner136_update) > 0
        sc = take!(chan)
        l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
        shell_edge = sc.shell_edge45

        updateShellEdge!(shell_edge, corner136_update[1], l₁)

        put!(chan, sc)
    end
    nothing
end

function updateShellEdge45!(chan::RemoteChannel{Channel{SubCuboid}}, edge36_update::Vector{Tuple{I,LatticeSite}},
                            corner236_update::Vector{LatticeSite}) where I<:Int
    
    sc = take!(chan)
    l₁ = sc.consts.l₁; l₂ = sc.consts.l₂; l₃ = sc.consts.l₃
    shell_edge = sc.shell_edge45

    updateShellEdge!(shell_edge, edge36_update)
    if length(corner236_update) > 0
        updateShellEdge!(shell_edge, corner236_update[1], 1)
    end

    put!(chan, sc)
    nothing
end

