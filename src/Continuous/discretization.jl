using IntervalMatrices: correction_hull

# ==================================
# Abstract interface
# ==================================

"""
    AbstractApproximationModel

Abstract supertype for all approximation models.
"""
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
_alias(setops::AbstractDirections) = setops
_alias(setops::Val{:lazy}) = setops
_alias(setops::Val{:concrete}) = setops
_alias(setops::Val{:vrep}) = setops

"""
    discretize(ivp::IVP, δ, alg::AbstractApproximationModel)

Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.

### Input

- `ivp`   -- initial value problem for a linear ODE in canonical form (see `Notes` below)
- `δ`     -- step size
- `alg`   -- algorithm used to compute the approximation model, choose among:

        - `Forward()`        -- use forward-time interpolation
        - `Backward()`      -- use backward-time interpolation
        - `CorrectionHull()` -- use the correction hull of the matrix exponential to a given order
        - `NoBloating()`     -- do not bloat the initial states

Different approximation algorithms and their respective options are described
in the docstring of each method, e.g. [`Forward`](@ref).

### Output

The initial value problem of a discrete system.

### Notes

Initial-value problems considered in this function are of the form

```math
x' = Ax(t) + u(t),\\qquad x(0) ∈ \\mathcal{X}_0,\\qquad (1)
```
and where ``u(t) ∈ U(k)`` add where ``\\{U(k)\\}_k`` is a sequence of sets of
non-deterministic inputs and ``\\mathcal{X}_0`` is the set of initial
states. Recall that this initial-value problem is called homogeneous whenever `U`
is the empty set. Other problems, e.g. ``x' = Ax(t) + Bu(t)`` can be brought
to the canonical form with the function [`normalize`](@ref).

The initial value problem returned by this function consists of a set discretized
(also called *bloated*) initial states ``Ω₀``, together with the coefficient matrix
``Φ = e^{Aδ}`` and a transformed sequence of inputs if ``U`` is non-empty.

Two main variations of this algorithm are considered: dense time case and
discrete time case.

- In the dense time case, the transformation is such that the trajectories
  of the given continuous system are included in the computed flowpipe of the
  discretized system. More precisely, given a step size ``δ`` and the system (1)
  conservative set-based discretization function computes a set, ``Ω₀``, that
  guarantees to contain all the trajectories of (1) starting at any ``x(0) ∈ \\mathcal{X}_0``
  and for any input function that satisfies ``u(t) ∈ U(1)``, for any ``t ∈ [0, δ]``.
  If ``U`` is time-varying, this function also discretizes the inputs for ``k ≥ 0``.

- In the discrete time case, there is no bloating of the initial states and the
  input is assumed to remain constant between sampled times. Use the algorithm
  `NoBloating()` for this setting. If ``U`` is time-varying, this function also discretizes
  the inputs for ``k ≥ 0``.

There are algorithms to obatin such transformations, called *approximation models*
in the technical literature. For references to the original papers, see the
docstring of each concrete subtype of `AbstractApproximationModel`.
"""
function discretize(ivp::IVP, δ, alg::AbstractApproximationModel) end

# ============================================================
# Approximation model in discrete time, i.e. without bloating
# ============================================================

"""
    NoBloating{EM, SO, IT} <: AbstractApproximationModel

No bloating, or discrete-time, approximation model.

### Fields

- `exp`     -- exponentiation method
- `setops`  -- set operations method
- `inv`     -- (optional, default: `false`) if `true`, assume that the state matrix
               is invertible and use its inverse in the `Φ` functions

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``
- ``Ω_0 ← \\mathcal{X}_0``
- ``V(k) ← Φ₁(A, δ)U(k)``, ``k ≥ 0``.

The function ``Φ₁(A, δ)`` is defined in [`Φ₁`](@ref).
We allow ``U`` to be a sequence of time varying non-deterministic input sets.

See also Eqs.(14) in [[BFFPSV18]](@ref).
"""
struct NoBloating{EM, SO, IT} <: AbstractApproximationModel
    exp::EM
    setops::SO
    inv::IT
end

# convenience constructor using symbols
function NoBloating(; exp=BaseExp, setops=:lazy, inv=false)
    return NoBloating(_alias(exp), _alias(setops), Val(inv))
end

function Base.show(io::IO, alg::NoBloating)
    print(io, "`NoBloating` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - invertibility assumption: $(alg.inv)")
end

Base.show(io::IO, m::MIME"text/plain", alg::NoBloating) = print(io, alg)

# homogeneous case
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::NoBloating)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)

    Ω0 = copy(X0) # alg.setops doesn't apply
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

# inhomogeneous case
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::NoBloating)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    U = next_set(inputset(ivp), 1)

    Φ = _exp(A, δ, alg.exp)
    if isa(U, AbstractSingleton)
        Mu = _Φ₁_u(A, δ, alg.exp, alg.inv, element(U), Φ)
        V = Singleton(Mu)
    else
        M = _Φ₁(A, δ, alg.exp, alg.inv, Φ)
        V = _apply_setops(M * U, alg.setops)
    end

    Ω0 = _initial_state(X0)

    In = IdentityMultiple(one(eltype(A)), size(A, 1))
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, In, X, V)
    return InitialValueProblem(Sdis, Ω0)
end

# extension for no bloating case and singleton initial states
_initial_state(X0::CartesianProduct{N, <:Singleton{N}, <:Singleton{N}}) where {N} = convert(Singleton, X0)
_initial_state(X0) = X0

# ==================================
# Forward approximation
# ==================================

"""
    Forward{EM, SO, SI, IT, BT} <: AbstractApproximationModel

Forward approximation model.

### Fields

- `exp`     -- exponentiation method
- `setops`  -- set opertaions method
- `sih`     -- symmetric interval hull
- `inv`     -- (optional, default: `false`) if `true`, assume that the state matrix
               is invertible and use its inverse in the `Φ` functions
- `backend` -- (optional, default: `nothing`) used if the algorithm needs to apply
               concrete polyhedral computations

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← CH(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU(0) ⊕ E_ψ(U(0), δ) ⊕ E^+(\\mathcal{X}_0, δ))``,
- ``V(k) ← δU(k) ⊕ E_ψ(U(k), δ)``.

Here we allow ``U`` to be a sequence of time varying non-deterministic input sets.

For the definition of the sets ``E_ψ`` and ``E^+`` see [[FRE11]](@ref).
The `Backward` method uses ``E^-``.
"""
struct Forward{EM, SO, SI, IT, BT} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
    inv::IT
    backend::BT
end

hasbackend(alg::Forward) = !isnothing(alg.backend)

# convenience constructor using symbols
function Forward(; exp=BaseExp, setops=:lazy, sih=:concrete, inv=false, backend=nothing)
    return Forward(_alias(exp), _alias(setops), Val(sih), Val(inv), backend)
end

function Base.show(io::IO, alg::Forward)
    print(io, "`Forward` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - symmetric interval hull method: $(alg.sih)\n")
    print(io, "    - invertibility assumption: $(alg.inv)")
    print(io, "    - polyhedral computations backend: $(alg.backend)")
end

Base.show(io::IO, m::MIME"text/plain", alg::Forward) = print(io, alg)

# ------------------------------------------------------------
# Forward Approximation: Homogeneous case
# ------------------------------------------------------------

# if A == |A|, then Φ can be reused in the computation of Φ₂(|A|, δ)
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = _Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)

    Einit = sih(P2A_abs * sih((A * A) * X0, alg.sih), alg.sih)
    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Einit)
    Ω0 = _apply_setops(Ω0, alg)
    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

function discretize(ivp::IVP{<:CLCS, Interval{N, IA.Interval{N}}}, δ, alg::Forward) where {N}
    A = state_matrix(ivp)
    @assert size(A, 1) == 1
    X0 = initial_state(ivp)

    a = A[1, 1]
    aδ = a * δ
    Φ = exp(aδ)
    A_abs = abs(a)

    # use inverse method
    @assert !iszero(a) "the given matrix should be invertible"

    # a_sqr = a * a
    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    #Einit = (P2A_abs * a_sqr) * RA._symmetric_interval_hull(X0).dat

    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    Einit = (Φ - one(N) - aδ) * _symmetric_interval_hull(X0).dat

    Ω0 = Interval(hull(X0.dat, Φ * X0.dat + Einit))
    X = stateset(ivp)
    # the system constructor creates a matrix
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# ------------------------------------------------------------
# Forward Approximation: Inhomogeneous case
# ------------------------------------------------------------

# TODO : specialize, add option to compute the concrete linear map
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = _Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)

    Einit = sih(P2A_abs * sih((A * A) * X0, alg.sih), alg.sih)

    U = next_set(inputset(ivp), 1)
    Eψ0 = sih(P2A_abs * sih(A * U, alg.sih), alg.sih)

    Ud = δ*U ⊕ Eψ0
    In = IdentityMultiple(one(eltype(A)), size(A, 1))

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Ud ⊕ Einit)
    Ω0 = _apply_setops(Ω0, alg.setops)
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, In, X, Ud)
    return InitialValueProblem(Sdis, Ω0)
end

# ==================================
# Backward approximation
# ==================================

"""
    Backward{EM, SO, SI, IT, BT} <: AbstractApproximationModel

Backward approximation model.

### Fields

- `exp`     -- exponentiation method
- `setops`  -- set opertaions method
- `sih`     -- symmetric interval hull
- `inv`     -- (optional, default: `false`) if `true`, assume that the state matrix
               is invertible and use its inverse in the `Φ` functions
- `backend` -- (optional, default: `nothing`) used if the algorithm needs to apply
               concrete polyhedral computations

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← CH(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU(0) ⊕ E_ψ(U(0), δ) ⊕ E^-(\\mathcal{X}_0, δ))``,
- ``V(k) ← δU(k) ⊕ E_ψ(U(k), δ)``.

Here we allow ``U`` to be a sequence of time varying non-deterministic input sets.

For the definition of the sets ``E_ψ`` and ``E^-`` see [[FRE11]](@ref).
The `Forward` method uses ``E^+``.
"""
struct Backward{EM, SO, SI, IT, BT} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
    inv::IT
    backend::BT
end

# convenience constructor using symbols
function Backward(; exp=BaseExp, setops=:lazy, sih=:concrete, inv=false, backend=nothing)
    return Backward(_alias(exp), _alias(setops), Val(sih), Val(inv), backend)
end

function Base.show(io::IO, alg::Backward)
    print(io, "`Backward` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - symmetric interval hull method: $(alg.sih)\n")
    print(io, "    - invertibility assumption: $(alg.inv)")
    print(io, "    - polyhedral computations backend: $(alg.backend)")
end

Base.show(io::IO, m::MIME"text/plain", alg::Backward) = print(io, alg)

# TODO: add corresponding `discrete` methods <<<<<

# ===============================================================
# Discretize using the correction hull of the matrix exponential
# ===============================================================

"""
    CorrectionHull{EM} <: AbstractApproximationModel

Discretization using the correction hull of the matrix exponential.

### Fields

- `exp`   -- exponentiation method
- `order` -- order of the Taylor series expansion of the matrix exponential

### Algorithm

For the homogeneous case, this method implements the transformation:

```math
Ω_0 = CH(X_0, e^{Aδ}  X_0) ⊕ FX_0
```
where ``F`` is the correction (interval) matrix.

For the inhomogeneous case, ``x' = Ax + u``,  ``x ∈ X, u ∈ U``, implements
``Ω_0 = CH(X_0, exp(Aδ)  X0) ⊕ FX0`` where ``F`` is the correction (interval) matrix.

In both cases, if ``A`` is an interval matrix, the exponential is overapproximated
using methods from `IntervalMatrices.jl`.
"""
struct CorrectionHull{EM} <: AbstractApproximationModel
   order::Int
   exp::EM
end

# convenience constructor using symbols
function CorrectionHull(; order::Int=10, exp=IntervalExpAlg(order))
    return CorrectionHull(order, _alias(exp))
end

function Base.show(io::IO, alg::CorrectionHull)
    print(io, "`CorrectionHull` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - order: $(alg.order)\n")
end

Base.show(io::IO, m::MIME"text/plain", alg::CorrectionHull) = print(io, alg)

# -----------------------------------------------------------------
# Correction hull: homogeneous case x' = Ax, x in X
# -----------------------------------------------------------------

function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::CorrectionHull)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)

    # compute exp(A*δ) * X0
    # TODO refactor / dispatch
    X0z = _convert_or_overapproximate(Zonotope, X0)

    Φ = _exp(A, δ, alg.exp)
    if isinterval(Φ)
        Y = _overapproximate(Φ * X0z, Zonotope)
    else
        Y = linear_map(Φ, X0z)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0 = minkowski_sum(H, R)

    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# -----------------------------------------------------------------
# Correction hull: inhomogeneous case x' = Ax + u, x in X, u ∈ U
# -----------------------------------------------------------------
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::CorrectionHull)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    U = next_set(inputset(ivp), 1) # inputset(ivp)
    n = size(A, 1)

    # here U is an interval matrix map of a lazyset, TODO refactor / dispatch
    if isa(U, LinearMap)
        Uz = _convert_or_overapproximate(Zonotope, LazySets.set(U))
        B = matrix(U)
        Uz = isinterval(B) ? _overapproximate(B * Uz, Zonotope) : linear_map(B, Uz)
    else # LazySet
        Uz = _convert_or_overapproximate(Zonotope, U)
    end
    if zeros(dim(U)) ∉ Uz
        error("this function is not implemented, see issue #253")
    end

    # TODO refactor Ω0_homog
    # TODO refactor / dispatch
    X0z = _convert_or_overapproximate(Zonotope, X0)

    Φ = _exp(A, δ, alg.exp)
    if isinterval(Φ)
        Y = _overapproximate(Φ * X0z, Zonotope)
    else
        Y = linear_map(Φ, X0z)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0_homog = minkowski_sum(H, R)

    # compute C(δ) * U
    Cδ = _Cδ(A, δ, alg.order)
    Ud = _overapproximate(Cδ * Uz, Zonotope)
    Ω0 = minkowski_sum(Ω0_homog, Ud)
    Idn = Φ # IntervalMatrix(one(A)) or IdentityMultiple(one(eltype(A)), n) # FIXME
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, Idn, X, Ud)
    return InitialValueProblem(Sdis, Ω0)
end

# =========================================================================
# Alternatives to apply the set operation depending on the desired output
# =========================================================================

function _apply_setops(X, alg::Forward)
    if hasbackend(alg)
        _apply_setops(X, alg.setops, alg.backend)
    else
        _apply_setops(X, alg.setops)
    end
end

_apply_setops(X::LazySet, ::Val{:lazy}) = X  # no-op
_apply_setops(X::LazySet, ::Val{:concrete}) = concretize(X)
_apply_setops(X, template::AbstractDirections) = overapproximate(X, template)
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
