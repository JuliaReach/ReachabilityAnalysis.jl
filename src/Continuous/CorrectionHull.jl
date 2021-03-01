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
    Φ = _exp(A, δ, alg.exp)

    # compute Ω0 = CH(X0, ΦX0) ⊕ FX0
    Ω0 = _discretize_chull(A, Φ, X0, δ, alg)

    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

function _discretize_chull(A, Φ::IntervalMatrix, X0, δ, alg)
    X0 = _convert_or_overapproximate(X0)
    Y = _overapproximate(Φ * X0, Zonotope)

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0 = minkowski_sum(H, R)
    return Ω0
end

# F(δ) without the E(δ) correction term
function _correction_hull_without_E(A, δ, p)
    timeint(δ, i) = interval((i^(-i / (i-1)) - i^(-1 / (i-1))) * δ^i, 0)
    F = sum(map(x -> timeint(δ, i) * x, A^i / factorial(i)) for i in 2:p)
    return IntervalMatrix(F)
end

function _discretize_chull(A, Φ::AbstractMatrix, X0, δ, alg)
    X0 = _convert_or_overapproximate(X0)
    Y = _overapproximate(Φ * X0, Zonotope)

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = _correction_hull_without_E(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0 = minkowski_sum(H, R)

#    Ω0_chull = minkowski_sum(convert(VPolygon, concretize(CH(X0, Φ*X0))),
#                  convert(VPolygon, FX0));

    return Ω0
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
