"""
    AbstractClusteringMethod{P}

Abstract supertype for all clustering types, with partition of type `P`.

### Notes

A clustering method defines a function which maps reach-sets to one or several
reach-sets in an over-approximative way, i.e. such that the set union of the input
reach-sets is included in the set union of the output reach-sets.

By taking the convex hull of the input reach-sets one can reduce the number of
outputs sets to a single one, overapproximately. This is the method that corresponds
to the `LazyClustering` type. However, in some cases it is conveinent to do other
types of transformations, such as:

- Return several reach-sets that are obtained by grouping the input in a way
  defined by a partition. For example, a vector of ten input sets can be clustered
  in two groups of five sets each, or five groups of two sets each, etc.
- Use different set representations such as boxes or zonotopes.
- Further post-process the output of the convexification by splitting into smaller
  sets, eg. box splitting.

Each concrete subtype of `AbstractClusteringMethod` has a parameter `P` that defines
what type of clustering strategy is applied. The method should be accessed with the
`partition` getter function.

The following strategies are implemented
at the interface level:

- If `P` is of type `Missing`: no partition is applied
- If `P` is of type integer: the partition corresponds to gruping the into the given integer
                             number of sets (or as close as possible)
- If `P` is of type vector of vectors: the given partition is applied
"""
abstract type AbstractClusteringMethod{P} end

_partition(C::AbstractClusteringMethod{Missing}, idx) = nfolds(idx, 1)
_partition(C::AbstractClusteringMethod{<:Integer}, idx) = nfolds(idx, partition(C))
_partition(C::AbstractClusteringMethod{VT}, idx) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = [idx[vi] for vi in partition(C)]

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

"""
    NoClustering{P} <: AbstractClusteringMethod{P}

No-op clustering method.
"""
struct NoClustering{P} <: AbstractClusteringMethod{P}
#
end

partition(::NoClustering{Missing}) = missing

NoClustering() = NoClustering(missing)

function cluster(F, idx, ::NoClustering)
    return view(F, idx)
end


# =====================================
# Lazy convexification clustering
# =====================================

"""
    LazyClustering{P} <: AbstractClusteringMethod{P}

Cluster according to the given partition by applying a lazy representation of the
convex hull.
"""
struct LazyClustering{P} <: AbstractClusteringMethod{P}
    partition::P
end

partition(C::LazyClustering) = C.partition

LazyClustering() = LazyClustering(missing)
LazyClustering(nchunks::D) where {D<:Integer} = LazyClustering{D}(nchunks)
LazyClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = LazyClustering{VT}(partition)

function cluster(F, idx, ::LazyClustering)
    return Convexify(view(F, idx))
end

# =====================================
# Box clustering
# =====================================

"""
    BoxClustering{PI, PO} <: AbstractClusteringMethod{P}

### Notes

This method first takes a lazy convex hull for the given partition, then computes
a tight hyperrectangular approximation for each element in the partition.
"""
struct BoxClustering{PI, PO} <: AbstractClusteringMethod{PI}
    inpartition::PI
    outpartition::PO
end

partition(C::BoxClustering) = C.inpartition

BoxClustering() = BoxClustering(missing, missing)
BoxClustering(nchunks::D) where {D<:Integer} = BoxClustering{D, Missing}(nchunks, missing)
BoxClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = BoxClustering{VT}(partition, missing)

# do not partition neither the input nor the output
function cluster(F::Flowpipe{N, ReachSet}, idx, ::BoxClustering{Missing, Missing}) where {N}
    convF = Convexify(view(F, idx))
    return [overapproximate(convF, Hyperrectangle)]
end

# partition the input array and do not partition the output
function cluster(F::Flowpipe{N, ReachSet}, idx, C::BoxClustering{PI, Missing}) where {N, PI}
    p = _partition(C, idx)
    convF = [Convexify(view(F, cj)) for cj in p]
    return [overapproximate(Xc, Hyperrectangle) for Xc in convF]
end

# we use a Zonotope overapproximation of the flowpipe, take their convex hull, and
# compute its box overapproximation
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N, ST}
    Fz = overapproximate(Zonotope, Flowpipe(view(F, idx)))
    return cluster(Fz, idx, C)
end

# TODO methods to split the output

#=
 # old
 function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N, ST}
     Fidx = Flowpipe(view(F, idx))
     charr = [set(overapproximate(Ri, Zonotope)) for Ri in Fidx]
     Xoa = overapproximate(ConvexHullArray(charr), Hyperrectangle)
     return [ReachSet(Xoa, tspan(Fidx))]
 end
=#

# =====================================
# Zonotope clustering
# =====================================

#=
struct ZonotopeClustering{P} <: AbstractClusteringMethod{P}
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

=#

# convexify and convert to vrep
#C = ReachabilityAnalysis.Convexify(sol[end-aux+1:end])
#Cvertex = convex_hull(vcat([vertices_list(Z) for Z in LazySets.array(set(C))]...)) |> VPolygon
