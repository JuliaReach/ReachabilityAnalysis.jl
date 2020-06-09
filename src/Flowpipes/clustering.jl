

# TODO document
abstract type AbstractClusteringMethod end


# partition an array x into n parts of equal size (if possible); if n is too large
# we split into length(x) parts
function nfolds(x::AbstractVector, n::Int)
    p = length(x)
    m = min(n, p)
    s = p / m
    [x[round(Int64, (i-1)*s)+1:min(p,round(Int64, i*s))] for i in 1:m]
end

# =====================================
# Identity (no clustering)
# =====================================

# TODO document
struct NoClustering <: AbstractClusteringMethod
#
end

# TODO document
function cluster(F, idx, ::NoClustering)
    return view(F, idx)
end

# =====================================
# Lazy convexification clustering
# =====================================

# TODO document
struct LazyClustering{P} <: AbstractClusteringMethod
    partition::P
end

LazyClustering() = LazyClustering(Val(false), missing)
LazyClustering(nchunks::Integer) = LazyClustering(Val(false), nchunks)
LazyClustering(partition::Vector{}) = LazyClustering(Val(true), partition)

function cluster(F, idx, ::LazyClustering)
    return Convexify(view(F, idx))
end


# =====================================
# Box clustering
# =====================================

"""
    BoxClustering{P} <: AbstractClusteringMethod

### Notes

This method first takes a lazy convex hull for the given partition, then computes
a tight hyperrectangular approximation for each element in the partition.
"""
struct BoxClustering{P} <: AbstractClusteringMethod
    partition::P
end

BoxClustering() = BoxClustering(missing)
BoxClustering(nchunks::D) where {D<:Integer} = BoxClustering{D}(nchunks)
BoxClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = BoxClustering{VT}(partition)

_partition(bc::BoxClustering{Missing}, idx) = nfolds(idx, 1)
_partition(bc::BoxClustering{<:Integer}, idx) = nfolds(idx, bc.partition)
_partition(bc::BoxClustering{VT}, idx) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = [idx[vi] for vi in bc.partition]

#= TEMP
# no partition
function cluster(F::Flowpipe{N, ReachSet}, idx, ::BoxClustering{Val{false}, Missing}) where {N}
    convF = Convexify(view(F, idx))
    return [overapproximate(convF, Hyperrectangle)]
end

# partition into the given number of chunks
function cluster(F::Flowpipe{N, ReachSet}, idx, bc::BoxClustering{Val{false}, <:Integer}) where {N}
    partition = nfolds(idx, bc.partition)
    convF = [Convexify(view(F, cj)) for cj in partition]
    return [overapproximate(Xc, Hyperrectangle) for Xc in convF]
end

partition(bc::Bo, idx)
=#

function cluster(F::Flowpipe{N, ReachSet}, idx, bc::BoxClustering{Val{true}}) where {N}
    convF = [Convexify(view(F, cj)) for cj in partition(bc, idx)]
    return [overapproximate(Xc, Hyperrectangle) for Xc in convF]
end

# we use a Zonotope overapproximation of the flowpipe, take thir convex hull, and
# compute its box overapproximation
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N, ST}
    Fidx = Flowpipe(view(F, idx))
    charr = [set(overapproximate(Ri, Zonotope)) for Ri in Fidx]
    Xoa = overapproximate(ConvexHullArray(charr), Hyperrectangle)
    return [ReachSet(Xoa, tspan(Fidx))]
end


# =====================================
# Zonotope clustering
# =====================================

struct ZonotopeClustering <: AbstractClusteringMethod
#
end

# for the generalization to > 2 ses, we iteratively apply the overapprox of the
# CH of two zonotopes, cf LazySets #2154
function cluster(F::Flowpipe{N, RT, VRT}, idx, ::ZonotopeClustering) where {N, ZT<:Zonotope, RT<:ReachSet{N, ZT}, VRT<:AbstractVector{RT}}
    if length(idx) == 1
        return F[idx]
    elseif length(idx) == 2
        X = ConvexHull(set(F[idx[1]]), set(F[idx[2]]))
        Y = overapproximate(X, Zonotope) # TODO pass algorithm
        return [ReachSet(Y, tspan(F[idx]))]
    else
        Zaux = overapproximate(ConvexHull(set(F[idx[1]]), set(F[idx[2]])), Zonotope)
        for k in 3:length(idx)
            Zaux = overapproximate(ConvexHull(Zaux, set(F[idx[k]])), Zonotope)
        end
        return [ReachSet(Zaux, tspan(F[idx]))]
    end
end

# convexify and convert to vrep
#C = ReachabilityAnalysis.Convexify(sol[end-aux+1:end])
#Cvertex = convex_hull(vcat([vertices_list(Z) for Z in LazySets.array(set(C))]...)) |> VPolygon
