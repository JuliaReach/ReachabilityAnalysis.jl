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

apply(f::DiscreteTransition,
      R::AbstractLazyReachSet,
      method::AbstractIntersectionMethod) = _apply(f, set(R), method)

# ==========================================
# Cases where the reset map is the identity
# ==========================================

# ---------------------------
# Intersection method: Exact
# ---------------------------

# all sets are polyhedral
function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::Exact) where {N,
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    clist_G = constraints_list(tr.G)
    clist_I⁻ = constraints_list(tr.I⁻)
    clist_I⁺ = constraints_list(tr.I⁺)
    clist_X = constraints_list(X)
    out = vcat(clist_G, clist_I⁻, clist_I⁺, clist_X)
    flag = remove_redundant_constraints!(out)
    if flag
        return HPolyhedron(out)
    else
        return EmptySet(dim(X))
    end
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::Exact) where {N,
                                GT<:AbstractPolyhedron{N},
                                VN,
                                IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                IT⁺<:AbstractPolyhedron{N}}

    clist_G = constraints_list(tr.G)
    clist_I⁺ = constraints_list(tr.I⁺)
    clist_X = constraints_list(X)
    aux = vcat(clist_G, clist_I⁺, clist_X)
    flag = remove_redundant_constraints!(aux)
    !flag && return EmptySet(dim(X))

    clist_I⁻ = tr.I⁻.array
    out = Vector{HPolyhedron{N, Vector{N}}}()

    for ci in clist_I⁻
        aux2 = vcat(aux, ci)
        flag = remove_redundant_constraints!(aux2)
        if !flag
            continue
        else
            push!(out, HPolyhedron(copy(aux2)))
        end
    end

    if isempty(out)
        return EmptySet(dim(X))
    else
        return UnionSetArray(out)
    end
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::ConvexHullArray{N, PT},
                ::Exact) where {N,
                                PT<:AbstractPolyhedron{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    vlist = [vertices_list(Xi) for Xi in array(X)]
    vlist_ch = convex_hull(vcat(vlist...)) # TODO pass the backend as an option
    Xhrep = tohrep(vlist_ch)
    out = _apply(tr, Xhrep, Exact())
    return HPolyhedron(out)
end

# ---------------------------------
# Intersection method: BoxEnclosure
# ---------------------------------

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::ConvexHullArray{N, PT},
                ::BoxEnclosure) where {N,
                                       PT<:AbstractPolyhedron{N},
                                       GT<:AbstractPolyhedron{N},
                                       VN,
                                       IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                       IT⁺<:AbstractPolyhedron{N}}
    Xbox = overapproximate(X, Hyperrectangle)
    out = _apply(tr, Xbox, Exact())
    return overapproximate(out, Hyperrectangle)
end

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::BoxEnclosure) where {N,
                                       GT<:AbstractPolyhedron{N},
                                       IT⁻<:AbstractPolyhedron{N},
                                       IT⁺<:AbstractPolyhedron{N}}

    out = _apply(tr, X, Exact())
    return overapproximate(out, Hyperrectangle)
end

function _apply(tr::DiscreteTransition{<:Universe, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolyhedron{N},
                ::BoxEnclosure) where {N,
                                       GT<:AbstractPolyhedron{N},
                                       VN,
                                       IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                       IT⁺<:AbstractPolyhedron{N}}
    out = _apply(tr, X, Exact())
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
    out = _apply(tr, Xbox, Exact())
    return overapproximate(out, Hyperrectangle)
end
