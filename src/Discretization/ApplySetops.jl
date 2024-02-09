# =========================================================================
# Alternatives to apply the set operation depending on the desired output
# =========================================================================
module ApplySetops

using ..DiscretizationModule: AbstractApproximationModel, hasbackend
using LazySets
using LazySets.Approximations: AbstractDirections

export _apply_setops

function _apply_setops(X, alg::AbstractApproximationModel)
    if hasbackend(alg)
        _apply_setops(X, alg.setops, alg.backend)
    else
        _apply_setops(X, alg.setops)
    end
end

_apply_setops(X::LazySet, ::Val{:lazy}) = X  # no-op
_apply_setops(X::LazySet, ::Val{:concrete}) = concretize(X) # concrete set
_apply_setops(X, template::AbstractDirections) = overapproximate(X, template) # template oa
_apply_setops(X::LazySet, ::Val{:box}) = box_approximation(X) # box oa

_apply_setops(M::AbstractMatrix, X::LazySet, ::Val{:lazy}) = M * X
_apply_setops(M::AbstractMatrix, X::LazySet, ::Val{:concrete}) = linear_map(M, X)

# evantually we should use concretize, but requires fast fallback operations in 2D
# such as Minkowski sum not yet available
function _apply_setops(X::ConvexHull{N,AT,MS}, ::Val{:vrep},
                       backend=nothing) where {N,
                                               AT<:AbstractPolytope{N},
                                               LM<:LinearMap{N,AT,N},
                                               MS<:MinkowskiSum{N,LM}}
    n = dim(X)
    VT = n == 2 ? VPolygon : VPolytope

    # CH(A, B) := CH(X₀, ΦX₀ ⊕ E₊)
    A = X.X
    B = X.Y
    X₀ = convert(VT, A)

    if n == 2
        ΦX₀ = convert(VT, B.X)
        E₊ = convert(VT, B.Y)
        out = convex_hull(X₀, minkowski_sum(ΦX₀, E₊))
    else
        # generic conversion to VPolytope is missing, see LazySets#2467
        ΦX₀ = VPolytope(vertices_list(B.X; prune=false))
        E₊ = convert(VT, B.Y)
        aux = minkowski_sum(ΦX₀, E₊; apply_convex_hull=false)
        out = convex_hull(X₀, aux; backend=backend)
    end

    return out
end

# give X = CH(X₀, ΦX₀ ⊕ E₊), return a zonotope overapproximation
function _apply_setops(X::ConvexHull{N,AT,MS}, ::Val{:zono},
                       backend=nothing) where {N,
                                               AT<:AbstractZonotope{N},
                                               LM<:LinearMap{N,AT,N},
                                               MS<:MinkowskiSum{N,LM}}
    # CH(A, B) := CH(X₀, ΦX₀ ⊕ E₊)
    A = X.X
    B = X.Y
    return overapproximate(CH(A, concretize(B)), Zonotope)
end

end  # module
