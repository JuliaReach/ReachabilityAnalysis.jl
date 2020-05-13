abstract type AbstractTransition end

"""
    DiscreteTransition{RT, GT, IT⁻, IT⁺, WT} <: AbstractTransition

Type that encodes a discrete transition with an affine assignment of the form:

```math
    post_d(X) = (R(X ∩ G ∩ I⁻) ⊕ W) ∩ I⁺
```
where ``I⁻`` and ``I⁺``  are invariants at the source and the target locations
respectively, ``G ⊆ \\mathbb{R}^n`` is the guard set of the transition, the assignment
is of the form ``x^+ := Rx + w``, ``w ∈ W``, ``x^+ ∈ \\mathbb{R}^m`` are the values after
the transition, ``R ∈ \\mathbb{R}^{m\\times n}`` is the assignment map (or reset map)
and ``W ⊆ \\mathbb{R}^m`` is a closed and bounded convex set of non-deterministic inputs.

### Fields

- `R`  -- assignment map
- `W`  -- non-deterministic inputs
- `G`  -- guard set of the transition from the source location to the target location
- `I⁻` -- invariant at the source location
- `I⁺` -- invariant at the target location
"""
struct DiscreteTransition{RT, WT, GT, IT⁻, IT⁺} <: AbstractTransition
    R::RT
    W::WT
    G::GT
    I⁻::IT⁻
    I⁺::IT⁺
end

reset_map(tr::DiscreteTransition{<:AbstractMatrix})  = ConstrainedAffineMap(tr.R, tr.W)
reset_map(tr::DiscreteTransition{<:IdentityMap}) = ConstrainedIdentityMap(dim(tr.W), tr.W)
guard(tr::DiscreteTransition) = tr.G
source_invariant(tr::DiscreteTransition) = tr.I⁻
target_invariant(tr::DiscreteTransition) = tr.I⁺

# constructor when the reset map is the identity
function DiscreteTransition(; guard::GT, source_invariant::IT⁻, target_invariant::IT⁺) where {GT, IT⁻, IT⁺}
    n = dim(source_invariant)
    return DiscreteTransition(IdentityMap(n), ZeroSet(n), guard, source_invariant, target_invariant)
end

# we assume that the reach-set is bounded: it's either a polytopic set
# or a union of polytopic sets
# => the return type is either the empty set or an HPolytope unless otherwise stated
# for examples in cases involving unions this method may return a union of polytopes
apply(tr::DiscreteTransition,
      R::AbstractLazyReachSet,
      method::AbstractIntersectionMethod) = _apply(tr, set(R), method)

# ==========================================
# Cases where the reset map is the identity
# ==========================================

# ----------------------------------------
# Intersection method: HRepIntersection
# ----------------------------------------

# all sets are polyhedral
function _apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolytope{N},
                method::HRepIntersection) where {N,
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    success, out = _intersection(X, tr.G, tr.I⁻, tr.I⁺, method)
    return success ? HPolytope(out) : EmptySet(dim(X))
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
function _apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolytope{N},
                method::HRepIntersection) where {N,
                                GT<:AbstractPolyhedron{N},
                                VN,
                                IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                IT⁺<:AbstractPolyhedron{N}}

    success, out = _intersection(X, tr.G, tr.I⁺, method)
    !success && return EmptySet(dim(X))

    clist_I⁻ = tr.I⁻.array
    out = Vector{HPolytope{N, Vector{N}}}()

    # intersect Y := X ∩ G ∩ I⁺ with each half-space ci in I⁻
    # only store those polytpes for which Y ∩ ci is non-empty
    for ci in clist_I⁻
        Y_ci = vcat(aux, ci)
        success = remove_redundant_constraints!(Y_ci)
        !success && continue
        push!(out, HPolytope(copy(Y_ci)))
    end

    return !isempty(out) ? UnionSetArray(out) : EmptySet(dim(X))
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
function _apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                X::ConvexHullArray{N, PT},
                method::HRepIntersection) where {N,
                                PT<:AbstractPolytope{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    vlist = vertices_list(X) # TODO pass the backend as an option
    Xhrep = tohrep(vlist)
    return  _apply(tr, Xhrep, method)
end

# ---------------------------------
# Intersection method: BoxEnclosure
# ---------------------------------
#
#= TODO: use SupportFunctionIntersection instead
#
function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::ConvexHullArray{N, PT},
                ::BoxEnclosure) where {N,
                                       PT<:AbstractPolyhedron{N},
                                       GT<:AbstractPolyhedron{N},
                                       VN,
                                       IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                       IT⁺<:AbstractPolyhedron{N}}
    Xbox = overapproximate(X, Hyperrectangle)
    out = _apply(tr, Xbox, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::BoxEnclosure) where {N,
                                       GT<:AbstractPolyhedron{N},
                                       IT⁻<:AbstractPolyhedron{N},
                                       IT⁺<:AbstractPolyhedron{N}}

    out = _apply(tr, X, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::BoxEnclosure) where {N,
                                       GT<:AbstractPolyhedron{N},
                                       VN,
                                       IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                       IT⁺<:AbstractPolyhedron{N}}
    out = _apply(tr, X, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::ConvexHullArray{N, PT},
                ::BoxEnclosure) where {N,
                                       PT<:AbstractPolyhedron{N},
                                       GT<:AbstractPolyhedron{N},
                                       IT⁻<:AbstractPolyhedron{N},
                                       IT⁺<:AbstractPolyhedron{N}}

    Xbox = overapproximate(X, Hyperrectangle)
    out = _apply(tr, Xbox, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

=#
