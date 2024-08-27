"""
    AbstractClusteringMethod{P}

Abstract supertype for all clustering types, with partition of type `P`.

### Notes

A clustering method defines a function which maps reach-sets to one or several
reach-sets. The mapping can possibly be over-approximative, i.e. such that the
set union of the input reach-sets is included in the set union of the output reach-sets.

By taking the convex hull of the input reach-sets one can reduce the number of
outputs sets to a single one, overapproximately. This is the method that corresponds
to the `LazyClustering` type. However, in some cases it is convenient to do other
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
- If `P` is of type integer: the partition corresponds to grouping the into the given integer
                             number of sets (or as close as possible)
- If `P` is of type vector of vectors: the given partition is applied

### Examples

`LazyClustering([1:2, 3:10])` groups the first two reach-sets in one cluster and
the third to tenth reach-sets in another cluster.
"""
abstract type AbstractClusteringMethod{P} end

_partition(::AbstractClusteringMethod{Missing}, idx) = nfolds(idx, 1)

_partition(method::AbstractClusteringMethod{Int}, idx) = nfolds(idx, partition(method))

function _partition(method::AbstractClusteringMethod{<:AbstractVector{<:AbstractVector{Int}}}, idx)
    return [idx[vi] for vi in partition(method)]
end

# partition an array x into n parts of equal size (if possible); if n is too large
# we split into length(x) parts
function nfolds(x::AbstractVector, n::Int)
    p = length(x)
    m = min(n, p)
    s = p / m
    return [x[(round(Int64, (i - 1) * s) + 1):min(p, round(Int64, i * s))] for i in 1:m]
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
struct LazyClustering{P,T} <: AbstractClusteringMethod{P}
    partition::P
    convex::T
end

partition(method::LazyClustering) = method.partition

LazyClustering(; convex::Bool=true) = LazyClustering(missing, Val(convex))
function LazyClustering(nchunks::Int; convex::Bool=true)
    return LazyClustering{Int,typeof(Val(convex))}(nchunks, Val(convex))
end
function LazyClustering(partition::VT) where {VT<:AbstractVector{<:AbstractVector{Int}}}
    return LazyClustering{VT,typeof(Val(convex))}(partition, Val(convex))
end

function cluster(F, idx, ::LazyClustering{Missing,Val{true}})
    return [convexify(view(F, idx))]
end

# TODO return a UnionSetArray of reach-sets
function cluster(F, idx, ::LazyClustering{Missing,Val{false}})
    return [view(F, idx)]
end

function cluster(F, idx, method::LazyClustering{P,Val{true}}) where {P}
    p = _partition(method, idx)
    return [convexify(view(F, cj)) for cj in p]
end

function cluster(F, idx, method::LazyClustering{P,Val{false}}) where {P}
    p = _partition(method, idx)
    return [view(F, cj) for cj in p]
end

# for Taylor model flowpipes we preprocess it with a zonotopic overapproximation
function cluster(F::Flowpipe{N,TaylorModelReachSet{N}}, idx,
                 method::LazyClustering{P,Val{true}}) where {N,P}
    return _cluster_TM(F, idx, method)
end
# disambiguations
function cluster(F::Flowpipe{N,TaylorModelReachSet{N}}, idx,
                 method::LazyClustering{Missing,Val{true}}) where {N}
    return _cluster_TM(F, idx, method)
end
function cluster(F::Flowpipe{N,TaylorModelReachSet{N}}, idx,
                 method::LazyClustering{P,Val{false}}) where {N,P}
    return _cluster_TM(F, idx, method)
end
function cluster(F::Flowpipe{N,TaylorModelReachSet{N}}, idx,
                 method::LazyClustering{Missing,Val{false}}) where {N}
    return _cluster_TM(F, idx, method)
end

function _cluster_TM(F, idx, method)
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)
    return cluster(Fz, 1:length(idx), method) # Fx is now indexed from 1 ... length(idx)
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
struct BoxClustering{PI,PO} <: AbstractClusteringMethod{PI}
    inpartition::PI
    outpartition::PO
end

partition(method::BoxClustering) = method.inpartition

BoxClustering() = BoxClustering(missing, missing)
BoxClustering(nchunks::Int) = BoxClustering{Int,Missing}(nchunks, missing)
function BoxClustering(partition::VT) where {VT<:AbstractVector{<:AbstractVector{Int}}}
    return BoxClustering{VT,Missing}(partition, missing)
end

_out_partition(::BoxClustering{PI,Missing}, idx, n) where {PI} = fill(1, n)
function _out_partition(method::BoxClustering{PI,Int}, idx, n) where {PI}
    return fill(method.outpartition, n)
end
function _out_partition(method::BoxClustering{PI,<:AbstractVector{Int}}, idx, n) where {PI}
    return method.outpartition
end

# do not partition neither the input nor the output
function cluster(F::Flowpipe{N,<:AbstractLazyReachSet}, idx,
                 ::BoxClustering{Missing,Missing}) where {N}
    convF = convexify(view(F, idx))
    return [overapproximate(convF, Hyperrectangle)]
end

# do not partition neither the input but split the output
function cluster(F::Flowpipe{N,<:AbstractLazyReachSet}, idx,
                 method::BoxClustering{Missing,PO}) where {N,PO}
    convF = convexify(view(F, idx))
    Y = overapproximate(convF, Hyperrectangle)
    p = _out_partition(method, idx, dim(F))
    return split(Y, p)
end

# partition the input array and do not partition the output
#function cluster(F::Flowpipe{N, Union{<:ReachSet, <:SparseReachSet}}, idx, method::BoxClustering{PI, Missing}) where {N, PI}
function cluster(F::Flowpipe{N,<:AbstractLazyReachSet}, idx,
                 method::BoxClustering{PI,Missing}) where {N,PI}
    p = _partition(method, idx)
    convF = [convexify(view(F, cj)) for cj in p]
    return [overapproximate(Xc, Hyperrectangle) for Xc in convF]
end

# partition the input array and do not partition the output
function cluster(F::Flowpipe{N,<:AbstractLazyReachSet}, idx,
                 method::BoxClustering{PI,PO}) where {N,PI,PO}
    p = _partition(method, idx)
    convF = [convexify(view(F, cj)) for cj in p]
    Y = [overapproximate(Xc, Hyperrectangle) for Xc in convF]
    pout = _out_partition(method, idx, dim(F))
    Ys = reduce(vcat, [split(y, pout) for y in Y])
    return Ys
end

# for Taylor model flowpipes we preprocess it with a zonotopic overapproximation
function cluster(F::Flowpipe{N,TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)

    # Fx is now indexed from 1 ... length(idx)
    return cluster(Fz, 1:length(idx), method)
end

# =====================================
# Box clustering with a vector output
# =====================================

struct BoxVecClustering{P} <: AbstractClusteringMethod{P}
    part::P
end

BoxVecClustering() = BoxVecClustering{Missing}(missing)

function cluster(F, idx, cl::BoxVecClustering)
    out = cluster(F, idx, BoxClustering(cl.part))
    return [out]
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
function cluster(F::Flowpipe{N,RT,VRT}, idx,
                 ::ZonotopeClustering{Missing}) where {N,ZT<:Zonotope,RT<:ReachSet{N,ZT},
                                                       VRT<:AbstractVector{RT}}
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

function _cartesian_product_standard(Z, blocks)
    @assert length.(blocks) == [2, 2, 1]

    # shuffle columns to bring to standard form
    M = Z.generators[:, [1, 2, 5, 6, 3, 4, 7, 8, 9]]

    # add extra column of zeros
    Mst = hcat(M[:, 1:4], zeros(5), M[:, 5:9])
    return Zst = Zonotope(Z.center, Mst)
end

# TODO wrap into clustering strategy and pass outer partition
function _cluster_decompose_zono(fp::AbstractVector{<:TaylorModelReachSet}, blocks, dirs)
    # zonotopic enclosure
    fpz = [overapproximate(R, Zonotope) for R in fp]
    return _cluster_decompose_zono(fpz, blocks, dirs)
end

function _cluster_decompose_zono(fpz::AbstractVector{RT}, blocks,
                                 dirs) where {N,ZT<:Zonotope,RT<:ReachSet{N,ZT}}

    # concrete projection onto each block
    πfpz = [[project(set(R), bi) for R in fpz] for bi in blocks]

    # reconstruct zonotope enclosure
    ZTD = Zonotope{Float64,Vector{Float64},Matrix{Float64}}
    zono_blk = Vector{ZTD}(undef, length(blocks))

    for (i, V) in enumerate(πfpz)
        Vch = ConvexHullArray(V)
        dim_blk = dim(Vch)
        @assert dim_blk == 1 || dim_blk == 2

        if dim_blk == 2
            # epsilon-close approximation
            Yi = overapproximate(Vch, 1e-3)

            # zonotope overapproximation
            Yi = overapproximate(Yi, Zonotope, dirs)

            # order reduction
            Yi = reduce_order(Yi, 2)

        elseif dim_blk == 1
            Yi = overapproximate(Vch, LazySets.Interval)
            Yi = convert(Zonotope, Yi)
        end

        Yi = Zonotope(Vector(Yi.center), Matrix(Yi.generators))
        zono_blk[i] = Yi
    end

    # recomposition into a higher-dimensional zonotope
    out = concretize(CartesianProductArray(zono_blk))

    # bring to standard form
    out_st = _cartesian_product_standard(out, blocks)

    return out_st
end

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
UnionClustering(nchunks::Int) = UnionClustering{Int}(nchunks)
function UnionClustering(partition::VT) where {VT<:AbstractVector{<:AbstractVector{Int}}}
    return UnionClustering{VT}(partition)
end

function cluster(F::Flowpipe{N,<:AbstractLazyReachSet}, idx, ::UnionClustering{Missing}) where {N}
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
function cluster(F::Flowpipe{N,<:TaylorModelReachSet{N}}, idx,
                 method::UnionClustering{P}) where {N,P}
    Fz = overapproximate(Flowpipe(view(F, idx)), Zonotope)

    # Fx is now indexed from 1 ... length(idx)
    return cluster(Fz, 1:length(idx), method)
end
