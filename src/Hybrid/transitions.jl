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
HybridSystems.guard(tr::DiscreteTransition) = tr.G
source_invariant(tr::DiscreteTransition) = tr.I⁻
target_invariant(tr::DiscreteTransition) = tr.I⁺

# constructor when the reset map is the identity
function DiscreteTransition(; guard::GT, source_invariant::IT⁻, target_invariant::IT⁺) where {GT, IT⁻, IT⁺}
    n = dim(source_invariant)
    return DiscreteTransition(IdentityMap(n), ZeroSet(n), guard, source_invariant, target_invariant)
end

# constructor from hybrid system given source and target modes id's and the transition t
function DiscreteTransition(H::HybridSystem, transition)
    assignment = resetmap(H, transition)
    R = _state_matrix(assignment)
    W = _affine_term(assignment)
    G = stateset(assignment) # HybridSystems.guard(H, transition) doesn't work
    I⁻ = stateset(H, source(H, transition))
    I⁺ = stateset(H, target(H, transition))
    return DiscreteTransition(R, W, G, I⁻, I⁺)
end

# -----------------------------------------------
# Getter functions for MathematicalSystems maps
# (see also MathematicalSystems#209)
# -----------------------------------------------

_state_matrix(m::IdentityMap) = m       # no-op, to dispatch on identity assignments
_state_matrix(m::ConstrainedIdentityMap) = IdentityMap(statedim(m))
_state_matrix(m::MathematicalSystems.LinearMap) = m.A
_state_matrix(m::ConstrainedLinearMap) = m.A
_state_matrix(m::MathematicalSystems.AffineMap) = m.A
_state_matrix(m::ConstrainedAffineMap) = m.A

_affine_term(m::IdentityMap) = ZeroSet(statedim(m)) # there is no numeric type in these maps; use ZeroSet default
_affine_term(m::ConstrainedIdentityMap) = ZeroSet(statedim(m))
_affine_term(m::MathematicalSystems.LinearMap{N}) where {N} = ZeroSet{N}(statedim(m))
_affine_term(m::ConstrainedLinearMap{N}) where {N} = ZeroSet{N}(statedim(m))
_affine_term(m::MathematicalSystems.AffineMap{N, <:AbstractVector}) where {N} = Singleton(m.c)
_affine_term(m::MathematicalSystems.AffineMap{N, <:LazySet}) where {N} = m.c
_affine_term(m::ConstrainedAffineMap{N, MT, VT, <:AbstractVector}) where {N, MT, VT} = Singleton(m.c)
_affine_term(m::ConstrainedAffineMap{N, MT, VT, <:LazySet}) where {N, MT, VT} = Singleton(m.c)

LazySets.dim(c::ConstrainedResetMap) = statedim(c)

function _state_matrix(m::Union{<:MathematicalSystems.ResetMap{N},
                                <:MathematicalSystems.ConstrainedResetMap{N}}) where {N}
    n = dim(m)
    v = ones(N, n)
    for i in keys(m.dict)
        v[i] = zero(N)
    end
    return Diagonal(v)
end

function _affine_term(m::Union{<:MathematicalSystems.ResetMap{N},
                               <:MathematicalSystems.ConstrainedResetMap{N}}) where {N}
    n = dim(m)
    b = sparsevec(Int[], N[], n)
    for (i, val) in m.dict
        b[i] = val
    end
    return Singleton(b)
end

# -------------------------------------------
# Applying reset maps
# -------------------------------------------

@inline function _apply_reset(R::AbstractMatrix, P::HPolytope, W::ZeroSet, ::HRepIntersection)
    return linear_map(R, P) # TODO pass options here
end

@inline function _apply_reset(R::AbstractMatrix, P::HPolytope, W::AbstractVector, ::HRepIntersection)
    return affine_map(R, P, W) # TODO pass options here
end

@inline function _apply_reset(R::AbstractMatrix, P::HPolytope, W::LazySet, ::HRepIntersection)
    Ylm = linear_map(R, P) # TODO pass options here
    return minkowski_sum(Ylm, W)
end

# we assume that the reach-set is bounded: it's either a polytopic set
# or a union of polytopic sets
# => the return type is either the empty set or an HPolytope unless otherwise stated
# for examples in cases involving unions this method may return a union of polytopes
apply(tr::DiscreteTransition,
      R::AbstractLazyReachSet,
      method::AbstractIntersectionMethod) = apply(tr, set(R), method)

# ==========================================
# Intersection method: HRepIntersection
# ==========================================

# -------------------------------------------
# Cases where the reset map is the identity
# -------------------------------------------

# all sets are polyhedral
function  apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolytope{N},
                method::HRepIntersection) where {N,
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    success, out = _intersection(X, tr.G, tr.I⁻, tr.I⁺, method)
    return success ? HPolytope(out) : EmptySet(dim(X))
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
               X::AbstractPolytope{N},
               method::HRepIntersection) where {N,
                                GT<:AbstractPolyhedron{N},
                                VN,
                                IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                IT⁺<:AbstractPolyhedron{N}}

    success, aux = _intersection(X, tr.G, tr.I⁺, method)
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
    return isempty(out) ? EmptySet(dim(X)) : UnionSetArray(out)
end

# the source invariant is the set union of half-spaces, the other sets are polyhedral
# and the set X is the convex hull of polytopes
function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
               X::ConvexHullArray{N, PT},
               method::HRepIntersection) where {N,
                                PT<:AbstractPolytope{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    vlist = vertices_list(X) # TODO pass the backend as an option
    Xhrep = tohrep(vlist)
    return apply(tr, Xhrep, method)
end

# the guard is the set union of half-spaces, the source invariant is polyhedral,
# and the target invariant is universal
function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                  X::AbstractPolytope{N},
                  method::HRepIntersection) where {N, VN,
                                                   GT<:UnionSetArray{N, HalfSpace{N, VN}},
                                                   IT⁻<:AbstractPolyhedron{N},
                                                   IT⁺<:Universe{N}}

    success, aux = _intersection(X, tr.I⁻, method)
    !success && return EmptySet(dim(X))

    clist_G = tr.G.array
    out = Vector{HPolytope{N, Vector{N}}}()

    # intersect Y := X ∩ I⁻ with each half-space ci in G
    # only store those polytpes for which Y ∩ ci is non-empty
    for ci in clist_G
        Y_ci = vcat(aux, ci)
        success = remove_redundant_constraints!(Y_ci)
        !success && continue
        push!(out, HPolytope(copy(Y_ci)))
    end

    return !isempty(out) ? UnionSetArray(out) : EmptySet(dim(X))
end


# -------------------------------------------
# Cases where the reset map is RX ⊕ W
# -------------------------------------------

function _apply_hrep(R, W, G, I⁻, I⁺, X, method::HRepIntersection)
    success, aux = _intersection(X, G, I⁻, method)
    !success && return EmptySet(dim(X))

    P = HPolytope(aux)
    Y = _apply_reset(R, P, W, method)

    success, out = _intersection(Y, I⁺, method)
    return success ? HPolytope(out) : EmptySet(dim(X))
end

# cases with W = 0
function apply(tr::DiscreteTransition{<:AbstractMatrix, <:ZeroSet, GT, IT⁻, IT⁺},
               X::AbstractPolytope{N},
               method::HRepIntersection) where {N, GT<:AbstractPolyhedron{N},
                                                   IT⁻<:AbstractPolyhedron{N},
                                                   IT⁺<:AbstractPolyhedron{N}}
    return _apply_hrep(tr.R, tr.W, tr.G, tr.I⁻, tr.I⁺, X, method)
end

# cases where the inputs are deterministic (W is a vector)
function apply(tr::DiscreteTransition{<:AbstractMatrix, <:AbstractVector, GT, IT⁻, IT⁺},
               X::AbstractPolytope{N},
               method::HRepIntersection) where {N, GT<:AbstractPolyhedron{N},
                                                   IT⁻<:AbstractPolyhedron{N},
                                                   IT⁺<:AbstractPolyhedron{N}}
    return _apply_hrep(tr.R, tr.W, tr.G, tr.I⁻, tr.I⁺, X, method)
end

# cases where the inputs are non-deterministic (W is a set)
function apply(tr::DiscreteTransition{<:AbstractMatrix, <:LazySet, GT, IT⁻, IT⁺},
               X::AbstractPolytope{N},
               method::HRepIntersection) where {N, GT<:AbstractPolyhedron{N},
                                                   IT⁻<:AbstractPolyhedron{N},
                                                   IT⁺<:AbstractPolyhedron{N}}
    return _apply_hrep(tr.R, tr.W, tr.G, tr.I⁻, tr.I⁺, X, method)
end

# ==========================================
# Intersection method: BoxIntersection
# ==========================================

function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                X::AbstractPolytope{N},
                ::BoxIntersection) where {N,
                                       GT<:AbstractPolyhedron{N},
                                       VN,
                                       IT⁻<:UnionSetArray{N, LinearConstraint{N, VN}},
                                       IT⁺<:AbstractPolyhedron{N}}

    out = apply(tr, X, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

# the guard is the set union of half-spaces, the source invariant is polyhedral,
# and the target invariant is the set union of half-spaces
function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
                  X::AbstractPolytope{N},
                  method::BoxIntersection) where {N, VN,
                                                  GT<:UnionSetArray{N, HalfSpace{N, VN}},
                                                  IT⁻<:AbstractPolyhedron{N},
                                                  IT⁺<:Universe{N}}

    out = apply(tr, X, HRepIntersection())
    return overapproximate(out, Hyperrectangle)
end

# all sets are polyhedral the reach-set is a set union
# each successor in the set union X is computed using HRep, then we overapproximate
# with a box
function apply(tr::DiscreteTransition{RT, <:ZeroSet, GT, IT⁻, IT⁺},
               X::UnionSetArray{N},
               method::BoxIntersection) where {N, RT,
                                                GT<:AbstractPolyhedron{N},
                                                IT⁻<:AbstractPolyhedron{N},
                                                IT⁺<:AbstractPolyhedron{N}}
    out = Vector{HPolytope{N, Vector{N}}}()
    for Xi in X.array
        # compute Yi := X ∩ G ∩ I⁻ exactly
        success, clist = _intersection(Xi, tr.G, tr.I⁻, HRepIntersection())
        !success && return EmptySet(dim(X))
        Yi = HPolytope(clist)

        # compute Ki := R * Yi exactly
        Ki = linear_map(tr.R, Yi)

        # compute Qi := Ki ∩ I⁺ exactly
        success, clist = _intersection(Ki, tr.I⁺, HRepIntersection())
        !success && return EmptySet(dim(X))

        # store the polytope in constraint representation
        Qi = HPolytope(clist)
        push!(out, Qi)
    end

    # overapproximate the result with a hyperrectangle
    return overapproximate(UnionSetArray(out), Hyperrectangle)
end

# all sets are polyhedral the reach-set is a set union
# each successor in the set union X is computed overapproximating with a
# hyperrectangle
#=
function apply(tr::DiscreteTransition{RT, <:ZeroSet, GT, IT⁻, IT⁺},
               X::UnionSetArray{N},
               method::BoxIntersection) where {N, RT,
                                                GT<:AbstractPolyhedron{N},
                                                IT⁻<:AbstractPolyhedron{N},
                                                IT⁺<:AbstractPolyhedron{N}}
    n = dim(X)
    boxdirs = BoxDirections(n)
    #out = Vector{Hyperrectangle{N, Vector{N}, Vector{N}}}()
    out = Vector{HPolytope{N, Vector{N}}}()
    for Xi in X.array
        # compute Yi := X ∩ G ∩ I⁻ using box directions
        Yi = _intersection(Xi, tr.G, tr.I⁻, method)

        # compute Ki := R * Yi using box directions
        Ki = box_approximation(tr.R * Yi)
        isempty(Ki) && continue

        # compute Qi := Ki ∩ I⁺ using box directions
        Qi = _intersection(Ki, tr.I⁺, method)

        # store the hyperrectangle
        push!(out, Qi)
    end

    # overapproximate the result with a hyperrectangle
    return overapproximate(UnionSetArray(out), Hyperrectangle)
end
=#

#=
function _apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
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
=#

#=


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

# ===============================================
# Intersection method: TemplateHullIntersection
# ===============================================

# compute (R(X ∩ G ∩ I⁻) ⊕ W) ∩ I⁺ first computing the exact intersection
# of Y := X ∩ G ∩ I⁻ all of which are polyhedral, using their list of constraints
# then we compute (RY ⊕ W) ∩ I⁺ using again the given template, Z := [(RY ⊕ W) ∩ I⁺]_dirs
function apply(tr::DiscreteTransition{<:AbstractMatrix, <:LazySet, GT, IT⁻, IT⁺},
               X::PT,
               method::TemplateHullIntersection) where {N,
                                PT<:AbstractPolyhedron{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    # intersect X with the guard and the source invariant (NOTE can be cached)
    success, Y = _intersection(X, tr.G, tr.I⁻, HRepIntersection())
    !success && return EmptySet(dim(X))

    # compute the intersection [(RY ⊕ W) ∩ I⁺]_dirs using the template
    K = (tr.R * HPolyhedron(Y)) ⊕ tr.W # lazy affine map

    Z = _intersection(K, tr.I⁺, method)
    return Z
end

# compute (R(X ∩ G ∩ I⁻) ⊕ W) ∩ I⁺ first computing the exact intersection
# of Y := G ∩ I⁻ all of which are polyhedral, then compute K := [X ∩ Y]_dirs,
# finally we compute (RK ⊕ W) ∩ I⁺ using again the given template, Z := [(RK ⊕ W) ∩ I⁺]_dirs
function apply(tr::DiscreteTransition{<:AbstractMatrix, <:LazySet, GT, IT⁻, IT⁺},
               X::PT,
               method::TemplateHullIntersection) where {N,
                                PT<:LazySet{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    # intersect the guard and the source invariant
    success, Y = _intersection(tr.G, tr.I⁻, HRepIntersection())
    !success && return EmptySet(dim(X))

    K = _intersection(X, HPolyhedron(Y), method)
    isempty(K) && return EmptySet(dim(X))

    # compute Z := [(RK ⊕ W) ∩ I⁺]_dirs using the template
    Km = (tr.R * K) ⊕ tr.W # lazy affine map
    Z = _intersection(Km, tr.I⁺, method)
    return Z
end

# compute X ∩ G ∩ I⁻ ∩ I⁺ first computing the exact intersection
# of Y := G ∩ I⁻ ∩ I⁺ all of which are polyhedral, then compute K := [X ∩ Y]_dirs
function apply(tr::DiscreteTransition{<:IdentityMap, <:ZeroSet, GT, IT⁻, IT⁺},
               X::PT,
               method::TemplateHullIntersection) where {N,
                                PT<:LazySet{N},
                                GT<:AbstractPolyhedron{N},
                                IT⁻<:AbstractPolyhedron{N},
                                IT⁺<:AbstractPolyhedron{N}}

    # intersect the guard and the source invariant
    success, Y = _intersection(tr.G, tr.I⁻, tr.I⁺, HRepIntersection())
    !success && return EmptySet(dim(X))

    # compute the intersection K := [X ∩ Y]_dirs using the template
    K = _intersection(X, HPolyhedron(Y), method)
    return K
end

# all sets are polyhedral the reach-set is a set union
# each successor in the set union X is computed using the template hull,
# then we overapproximate the union again using the template hull
function apply(tr::DiscreteTransition{RT, WT, GT, IT⁻, IT⁺},
               X::UnionSetArray{N},
               method::TemplateHullIntersection{N, VN, TN}) where {N, RT<:AbstractMatrix{N}, WT<:LazySet{N},
                                                        GT<:AbstractPolyhedron{N},
                                                        IT⁻<:AbstractPolyhedron{N},
                                                        IT⁺<:AbstractPolyhedron{N}, VN, TN}

    # compute Y := X ∩ G ∩ I⁻ using the template
    Y = _intersection(X, tr.G, tr.I⁻, method)

    # compute K := (R*Y⊕W) ∩ I⁺ using the template
    K = _intersection(tr.R * Y ⊕ tr.W, tr.I⁺, method)

    return K
end
