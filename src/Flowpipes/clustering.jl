"""
    AbstractClusteringMethod{P}

Abstract supertype for all clustering types, with partition of type `P`.

### Notes

A clustering method defines a function which maps reach-sets to one or several
reach-sets. The mapping can possibly be over-approximative, i.e. such that the
set union of the input reach-sets is included in the set union of the output reach-sets.

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

### Examples

`LazyClustering([1:2, 3:10])` groups the first two reach-sets in one cluster and
the third to tenth reach-sets in another cluster.
"""
abstract type AbstractClusteringMethod{P} end

_partition(method::AbstractClusteringMethod{Missing}, idx) = nfolds(idx, 1)
_partition(method::AbstractClusteringMethod{<:Integer}, idx) = nfolds(idx, partition(method))
_partition(method::AbstractClusteringMethod{VT}, idx) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = [idx[vi] for vi in partition(method)]

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

NoClustering() = NoClustering{Missing}()

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
struct LazyClustering{P, T} <: AbstractClusteringMethod{P}
    partition::P
    convex::T
end

partition(method::LazyClustering) = method.partition

LazyClustering(; convex::Bool=true) = LazyClustering(missing, Val(convex))
LazyClustering(nchunks::D; convex::Bool=true) where {D<:Integer} = LazyClustering{D, typeof(Val(convex))}(nchunks, Val(convex))
LazyClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = LazyClustering{VT, typeof(Val(convex))}(partition, Val(convex))

function cluster(F, idx, ::LazyClustering{Missing, Val{true}})
    return [convexify(view(F, idx))]
end

# FIXME return a UnionSetArray of reach-sets
function cluster(F, idx, ::LazyClustering{Missing, Val{false}})
    return [view(F, idx)]
end

function cluster(F, idx, method::LazyClustering{P, Val{true}}) where {P}
    p = _partition(method, idx)
    convF = [convexify(view(F, cj)) for cj in p]
end

function cluster(F, idx, method::LazyClustering{P, Val{false}}) where {P}
    p = _partition(method, idx)
    convF = [view(F, cj) for cj in p]
end

# for Taylor model flowpipes we preprocess it with a zonotopic overapproximation
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::LazyClustering{P, Val{true}}) where {N, P}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)
    return cluster(Fz, 1:length(idx), method) # Fx is now indexed from 1 ... length(idx)
end

# ambiguity fix
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::LazyClustering{P, Val{false}}) where {N, P}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)
    return cluster(Fz, 1:length(idx), method)
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

partition(method::BoxClustering) = method.inpartition

BoxClustering() = BoxClustering(missing, missing)
BoxClustering(nchunks::D) where {D<:Integer} = BoxClustering{D, Missing}(nchunks, missing)
BoxClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = BoxClustering{VT, Missing}(partition, missing)

_out_partition(method::BoxClustering{PI, Missing}, idx, n) where {PI} = fill(1, n)
_out_partition(method::BoxClustering{PI, <:Integer}, idx, n) where {PI} = fill(method.outpartition, n)
_out_partition(method::BoxClustering{PI, VT}, idx, n) where {PI, D<:Integer, VT<:AbstractVector{D}} = method.outpartition

# do not partition neither the input nor the output
function cluster(F::Flowpipe{N, <:AbstractLazyReachSet}, idx, ::BoxClustering{Missing, Missing}) where {N}
    convF = convexify(view(F, idx))
    return [overapproximate(convF, Hyperrectangle)]
end

# do not partition neither the input but split the output
function cluster(F::Flowpipe{N, <:AbstractLazyReachSet}, idx, method::BoxClustering{Missing, PO}) where {N, PO}
    convF = convexify(view(F, idx))
    Y = overapproximate(convF, Hyperrectangle)
    p = _out_partition(method, idx, dim(F))
    return split(Y, p)
end

# partition the input array and do not partition the output
#function cluster(F::Flowpipe{N, Union{<:ReachSet, <:SparseReachSet}}, idx, method::BoxClustering{PI, Missing}) where {N, PI}
function cluster(F::Flowpipe{N, <:AbstractLazyReachSet}, idx, method::BoxClustering{PI, Missing}) where {N, PI}
    p = _partition(method, idx)
    convF = [convexify(view(F, cj)) for cj in p]
    return [overapproximate(Xc, Hyperrectangle) for Xc in convF]
end

# partition the input array and do not partition the output
function cluster(F::Flowpipe{N, <:AbstractLazyReachSet}, idx, method::BoxClustering{PI, PO}) where {N, PI, PO}
    p = _partition(method, idx)
    convF = [convexify(view(F, cj)) for cj in p]
    Y = [overapproximate(Xc, Hyperrectangle) for Xc in convF]
    pout = _out_partition(method, idx, dim(F))
    Ys = reduce(vcat, [split(y, pout) for y in Y])
    return Ys
end

# for Taylor model flowpipes we preprocess it with a zonotopic overapproximation
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)

    # Fx is now indexed from 1 ... length(idx)
    return cluster(Fz, 1:length(idx), method)
end

# =====================================
# Zonotope clustering
# =====================================

"""
    ZonotopeClustering{P} <: AbstractClusteringMethod{P}

### Notes

This method first takes a lazy convex hull for the given partition, then computes
a zonotope overapproximation of the convex hull.
"""
struct ZonotopeClustering{P} <: AbstractClusteringMethod{P}
    partition::P
end

ZonotopeClustering() = ZonotopeClustering(missing)

# for the generalization to > 2 ses, we iteratively apply the overapprox of the
# CH of two zonotopes, cf LazySets #2154
function cluster(F::Flowpipe{N, RT, VRT}, idx, ::ZonotopeClustering{Missing}) where {N, ZT<:Zonotope, RT<:ReachSet{N, ZT}, VRT<:AbstractVector{RT}}
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
#C = ReachabilityAnalysis.convexify(sol[end-aux+1:end])
#Cvertex = convex_hull(vcat([vertices_list(Z) for Z in LazySets.array(set(C))]...)) |> VPolygon

# =====================================
# Lazy union set array clustering
# =====================================

"""
    UnionClustering{P} <: AbstractClusteringMethod{P}

Cluster according to the given partition by applying a lazy representation of the
set union.
"""
struct UnionClustering{P} <: AbstractClusteringMethod{P}
    partition::P
end

partition(method::UnionClustering) = method.partition

UnionClustering() = UnionClustering(missing)
UnionClustering(nchunks::D) where {D<:Integer} = UnionClustering{D}(nchunks)
UnionClustering(partition::VT) where {D<:Integer, VTi<:AbstractVector{D}, VT<:AbstractVector{VTi}} = UnionClustering{VT}(partition)

function cluster(F::Flowpipe{N, <:AbstractLazyReachSet}, idx, method::UnionClustering{Missing}) where {N}
    Fidx = view(F, idx)
    Δt = tspan(Fidx)
    Uidx = UnionSetArray([set(R) for R in Fidx])
    return [ReachSet(Uidx, Δt)]
end

function cluster(F::Flowpipe, idx, method::UnionClustering{P}) where {P}
    p = _partition(method, idx)
    return [cluster(F, cj, UnionClustering()) for cj in p]
end

# for Taylor model flowpipes we preprocess it with a zonotopic overapproximation
function cluster(F::Flowpipe{N, <:TaylorModelReachSet{N}}, idx, method::UnionClustering{P}) where {N, P}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)

    # Fx is now indexed from 1 ... length(idx)
    return cluster(Fz, 1:length(idx), method)
end
