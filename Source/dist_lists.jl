####################################################################################################
#                            Implementation of Distributed Lists
#
####################################################################################################

# The list gets divided into chuncks. Each Chunck has a p which stores the pid of the process where
# the chunck is stored, a chan that is a remote handle to the array chunck stored in the chunck
# and the part of the original range that is stored in the chunck.
@everywhere struct DLChunck{T}
    p::Int64                    
    chan::RemoteChannel{Channel{Array{T,1}}}
    range::UnitRange{Int64}
end

# Each DList then essentially is an array of chuncks.
struct DList
    chunck_list::Array{DLChunck, 1}
    elem_pr_chunck::Array{Int64, 1}               # Number of elements stored in chunck i.
    chuncks::Int64                                # Number of chuncks in chunck_list
    chunck_map::Array{Tuple{Int64, Int64}, 1}     # Map from the original indices to a tuple consisting of the index
    # of the chunck containing the index as well as the position in the chunck.
end

# Hack for fixing remote channels in 1.0.3 commit 099e826241
@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()

# --------------------------------------------------------------------------------------------------
import Base.length
function length(dlist::DList)
    l = 0
    for i = 1:dlist.chuncks
        l += dlist.elem_pr_chunck[i]
    end
    return l
end        

# --------------------------------------------------------------------------------------------------
# Divides the number of elements in the list n into evenly divided chunck-ranges over the available workers.
function getChunckRanges(n::T, nw::Int64=nprocs()-1) where T <: Int
    chunck_min, pluss_num = divrem(n, nw)
    chunck_ranges = Array{UnitRange{Int64}, 1}(undef, nw)
    for i = 1:pluss_num
        chunck_ranges[i] = range((i-1)*(chunck_min+1)+1; length=chunck_min+1)
    end
    for i = pluss_num+1:nw
        chunck_ranges[i] = range(pluss_num*(chunck_min+1)+(i-1-pluss_num)*chunck_min+1; length=chunck_min)
    end
    return chunck_ranges
end

# --------------------------------------------------------------------------------------------------
# Uses the chunck-ranges to copy the list-elements to the different processes for each chunck.
function distribute(r_list::Array{T,1}; pids = workers()) where {T}
    n = length(r_list)
    nw = length(pids)
    chunck_ranges = getChunckRanges(n, nw)
    
    chunck_list = Array{DLChunck, 1}(undef, nw)
    for (i, p) = enumerate(pids)
        chunck_list[i] = DLChunck(p, RemoteChannel(()->Channel{Array{T,1}}(1), p), chunck_ranges[i])
        put!(chunck_list[i].chan, r_list[chunck_ranges[i]])
    end
    chunck_list
end

# --------------------------------------------------------------------------------------------------
# Constructor for DList. Takes an original list r_list and divides it into chuncks over the available
# workers.
function DList(r_list::Array{T, 1}) where {T}
    chunck_list = distribute(r_list)
    chuncks = length(chunck_list)
    elem_pr_chunck = Array{Int64, 1}(undef, chuncks)
    chunck_map = Array{Tuple{Int64, Int64}, 1}(undef, length(r_list))
    k = 1
    for i = 1:chuncks
        elem_pr_chunck[i] = length(fetch(chunck_list[i].chan))
        for j = 1:elem_pr_chunck[i]
            chunck_map[k] = (i, j)
            k += 1
        end
    end
    DList(chunck_list, elem_pr_chunck, chuncks, chunck_map)
end

# --------------------------------------------------------------------------------------------------
# Takes a function f which is assumed to take an element of the original list as only argument, applies
# this function to all element in list and output an array of the results.
@everywhere function dmap(f::Function, chan::RemoteChannel{Channel{Array{T, 1}}}) where {T}
    local_list = fetch(chan)
    [f(el) for el in local_list]
end
function dmap(f::Function, dlist::DList)
    chuncks = dlist.chuncks
    futures = Array{Future, 1}(undef, chuncks)
    
    for (i, ck) = enumerate(dlist.chunck_list)
        futures[i] = @spawnat ck.p dmap(f, ck.chan)
    end
    vcat([fetch(futures[i]) for i = 1:chuncks]...)
end

# --------------------------------------------------------------------------------------------------
# Takes a function f! assumed to take an element of the original list as only argument and applies
# this function to all elements of the list.
@everywhere function dmutate(f!::Function, channel::RemoteChannel{Channel{Array{T, 1}}}) where {T}
    local_list = take!(channel)
    for el in local_list
        f!(el)
    end
    put!(channel, local_list)
    nothing
end
function dmutate(f!::Function, list::DList)
    chuncks = list.chuncks
    futures = Array{Future, 1}(undef, chuncks)
    
    for (i, ck) = enumerate(list.chunck_list)
        futures[i] = @spawnat ck.p dmutate(f!, ck.chan)
    end
    for i = 1:chuncks
        wait(futures[i])
    end
    nothing
end

# --------------------------------------------------------------------------------------------------
# Same as dmap, but now the array elements can be changed by f!.
@everywhere function impureMutate(f!::Function, chan::RemoteChannel{Channel{Array{T, 1}}}) where {T}
    local_list = take!(chan)
    res_list = [f!(el) for el in local_list]
    put!(chan, local_list)
    res_list
end
function impureMutate(f!::Function, dlist::DList)
    chuncks = dlist = chuncks
    futures = Array{Future, 1}(undef, chuncks)
    
    for (i, ck) = enumerate(dlist.chunck_list)
        futures[i] = @spawnat ck.p impureMutat(f!, ck.chan)
    end
    vcat([fetch(futures[i]) for i = 1:chuncks]...)
end

# --------------------------------------------------------------------------------------------------
# Brings a copy of the whole list into the calling process.
function localize(chunck_list::Array{DLChunck, 1})
    return vcat([fetch(ck.chan) for ck in chunck_list]...)
end
function localize(dlist::DList)
    localize(dlist.chunck_list)
end
