using IntervalMatrices: correction_hull, input_correction

# ==================================
# Abstract interface
# ==================================

abstract type AbstractApproximationModel end

function _default_approximation_model(ivp::IVP{<:AbstractContinuousSystem})
    return Forward()
end

# some algorithms require a polyhedral computations backend
hasbackend(alg::AbstractApproximationModel) = false

# symmetric inteval hull options
sih(X, ::Val{:lazy}) = SymmetricIntervalHull(X)
sih(X, ::Val{:concrete}) = _symmetric_interval_hull(X)

# interval matrix functions
isinterval(A::AbstractMatrix{N}) where {N<:Number} = false
isinterval(A::IntervalMatrix{N, IT}) where {N, IT<:IA.Interval{N}} = true
isinterval(A::AbstractMatrix{IT}) where {IT<:IA.Interval} = true

# options for a-posteriori transformation of a discretized set
# valid options are:
# AbstractDirections, Val{:lazy}, Val{:concrete}, Val{:vrep}, Val{:zono}, Val{:zonotope}
# _alias(setops) = setops # no-op

_alias(setops::AbstractDirections) = setops
_alias(setops::Val{:lazy}) = setops
_alias(setops::Val{:concrete}) = setops
_alias(setops::Val{:vrep}) = setops
_alias(setops::Val{:box}) = setops
_alias(setops::Val{:zono}) = setops
_alias(setops::Val{:zonotope}) = Val(:zono)

"""
    discretize(ivp::IVP, δ, alg::AbstractApproximationModel)

Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.

### Input

- `ivp`   -- initial value problem for a linear ODE in canonical form (see `Notes` below)
- `δ`     -- step size
- `alg`   -- algorithm used to compute the approximation model

### Output

The initial value problem of a discrete system.

### Notes

Different approximation algorithms and their respective options are described
in the docstring of each method, e.g. [`Forward`](@ref). Here is a list of all
the available approximation models:

```jldoctest
julia> subtypes(ReachabilityAnalysis.AbstractApproximationModel)
5-element Vector{Any}:
 Backward
 CorrectionHull
 Forward
 NoBloating
 StepIntersect
```

Initial-value problems considered in this function are of the form

```math
x' = Ax(t) + u(t),\\qquad x(0) ∈ \\mathcal{X}_0,\\qquad (1)
```
and where ``u(t) ∈ U(k)`` add where ``\\{U(k)\\}_k`` is a sequence of sets of
non-deterministic inputs and ``\\mathcal{X}_0`` is the set of initial
states. Other problems, e.g. ``x' = Ax(t) + Bu(t)`` can be brought
to the canonical form with the function [`normalize`](@ref).

For references to the original papers introducing each algorithm, see the docstrings,
e.g. `?Forward`.
"""
function discretize(ivp::IVP, δ, alg::AbstractApproximationModel)
    error("discretization not implemented for the given arguments: $ivp, $alg")
end

# =========================================
# Conservative time discretization methods
# =========================================

# Approximation model in discrete time, i.e. without bloating
include("NoBloating.jl")

# Forward approximation
include("Forward.jl")

# Backward approximation
include("Backward.jl")

# Intersect one step forward in time with one step backward
include("StepIntersect.jl")

# Discretize using the correction hull of the matrix exponential
include("CorrectionHull.jl")

# First-order approximation from d/dt
include("SecondOrderddt.jl")

# First-order approximation with zonotope
include("FirstOrderZonotope.jl")

# First-order approximation
include("FirstOrder.jl")

# Forward-Backward discretization using continuous convex hull
include("ForwardBackward.jl")

# =========================================================================
# Alternatives to apply the set operation depending on the desired output
# =========================================================================

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
function _apply_setops(X::ConvexHull{N, AT, MS}, ::Val{:vrep}, backend=nothing) where {N,
                            AT<:AbstractPolytope{N}, MT,
                            LM<:LinearMap{N, AT, N, MT},
                            BT<:AbstractPolytope, MS<:MinkowskiSum{N, LM}}
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
        ΦX₀ = VPolytope(vertices_list(B.X, prune=false))
        E₊ = convert(VT, B.Y)
        aux = minkowski_sum(ΦX₀, E₊, apply_convex_hull=false)
        out = convex_hull(X₀, aux, backend=backend)
    end

    return out
end

# give X = CH(X₀, ΦX₀ ⊕ E₊), return a zonotope overapproximation
function _apply_setops(X::ConvexHull{N, AT, MS}, ::Val{:zono}, backend=nothing) where {N,
                            AT<:AbstractZonotope{N}, MT,
                            LM<:LinearMap{N, AT, N, MT},
                            BT<:AbstractZonotope, MS<:MinkowskiSum{N, LM}}
    # CH(A, B) := CH(X₀, ΦX₀ ⊕ E₊)
    A = X.X
    B = X.Y
    return overapproximate(CH(A, concretize(B)), Zonotope)
end
