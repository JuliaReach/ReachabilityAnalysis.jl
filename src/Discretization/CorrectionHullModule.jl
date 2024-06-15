# ===============================================================
# Discretize using the correction hull of the matrix exponential
# ===============================================================
module CorrectionHullModule

using ..DiscretizationModule
using ..Exponentiation: _exp, _alias, IntervalExpAlg
using ..Overapproximate: _convert_or_overapproximate, _overapproximate
using IntervalMatrices: AbstractIntervalMatrix, IntervalMatrix, correction_hull,
                        input_correction, _exp_remainder
using MathematicalSystems
using LazySets
using LazySets: LinearMap
import IntervalArithmetic as IA
using LinearAlgebra
using Reexport

export CorrectionHull

@reexport import ..DiscretizationModule: discretize

const CLCS = ConstrainedLinearContinuousSystem
const CLCCS = ConstrainedLinearControlContinuousSystem

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
    print(io, "`CorrectionHull` approximation model with:\n")
    print(io, "    - exponentiation method: $(alg.exp)\n")
    print(io, "    - order: $(alg.order)\n")
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", alg::CorrectionHull) = print(io, alg)

# -----------------------------------------------------------------
# Correction hull: homogeneous case x' = Ax, x in X
# -----------------------------------------------------------------

function discretize(ivp::IVP{<:CLCS,<:LazySet}, δ, alg::CorrectionHull)
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

function _discretize_chull(A, Φ::IntervalMatrix, X0, δ, alg, P=nothing)
    X0z = _convert_or_overapproximate(X0, Zonotope)
    Y = _overapproximate(Φ * X0z, Zonotope)
    if !isnothing(P)
        Y = minkowski_sum(Y, P)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0 = minkowski_sum(H, R)
    return Ω0
end

# F(δ) without the E(δ) correction term
function _correction_hull_without_E(A, δ, p)
    timeint(δ, i) = IA.interval((i^(-i / (i - 1)) - i^(-1 / (i - 1))) * δ^i, 0)
    F = sum(map(x -> timeint(δ, i) * x, A^i / factorial(i)) for i in 2:p)
    return IntervalMatrix(F)
end

function _discretize_chull(A, Φ::AbstractMatrix, X0, δ, alg)
    X0z = _convert_or_overapproximate(X0, Zonotope)
    Y = linear_map(Φ, X0z)

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = _correction_hull_without_E(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)

    Ω0 = minkowski_sum(H, R)
    return Ω0
end

# -----------------------------------------------------------------
# Correction hull: inhomogeneous case x' = Ax + u, x in X, u ∈ U
# -----------------------------------------------------------------
function discretize(ivp::IVP{<:CLCCS,<:LazySet}, δ, alg::CorrectionHull)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    U = next_set(inputset(ivp), 1) # inputset(ivp)

    # here U is an interval matrix map of a lazyset, TODO refactor / dispatch
    if isa(U, LinearMap)
        Uz = _convert_or_overapproximate(Zonotope, LazySets.set(U))
        B = matrix(U) # TODO remove this case, since preprocessing normalize(ivp) would make B the identity
        Uz = isinterval(B) ? _overapproximate(B * Uz, Zonotope) : linear_map(B, Uz)
    else # LazySet
        Uz = _convert_or_overapproximate(Zonotope, U)
    end

    Φ = _exp(A, δ, alg.exp)

    A_interval = A isa AbstractIntervalMatrix ? A : IntervalMatrix(A)

    origin_not_contained_in_U = zeros(dim(U)) ∉ Uz

    if origin_not_contained_in_U
        # shift U to origin
        u = center(Uz)
        Uz = Zonotope(zeros(dim(U)), genmat(Uz))
    end

    # note: this evaluates Φ₁(A, δ) using Taylor expansion of the interval matrix
    Cδ = _Cδ(A_interval, δ, alg.order)

    if origin_not_contained_in_U
        # compute C(δ) * u
        Pu = _convert_intervals_to_box_with_vectors(Cδ * u)
    else
        Pu = nothing
    end

    Ω0 = _discretize_chull(A_interval, Φ, X0, δ, alg, Pu)

    if origin_not_contained_in_U
        # apply another correction hull for the shifted U
        F = input_correction(A_interval, δ, alg.order)
        Fu = _convert_intervals_to_box_with_vectors(F * u)
        Ω0 = minkowski_sum(Ω0, Fu)
    end
    # Ω0 = _apply_setops(Ω0, alg.setops) # TODO requires to add `setops` field to the struct

    # compute C(δ) * U
    Ud = _overapproximate(Cδ * Uz, Zonotope)

    # add inputs to Ω0
    Ω0 = minkowski_sum(Ω0, Ud)

    B = Matrix(IdentityMultiple(one(eltype(A)), size(A, 1)))
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, B, X, Ud)
    return InitialValueProblem(Sdis, Ω0)
end

# convert a vector of intervals to a hyperrectangle with normal vector type
function _convert_intervals_to_box_with_vectors(intervals)
    # convert to Hyperrectangle
    H = convert(Hyperrectangle, IntervalBox(intervals))
    # remove static vectors (they scale poorly)
    H = Hyperrectangle(Vector(center(H)), Vector(radius_hyperrectangle(H)))
    return H
end

# TODO: outsource to IntervalMatrices.jl
# compute Iδ + 1/2 * δ^2 * A + 1/6 * δ^3 * A² + ... + 1/(η+1)! * δ^(η+1) * A^η + E(δ) * δ
function _Cδ(A, δ, order)
    n = size(A, 1)
    A² = A * A
    if isa(A, IntervalMatrix)
        Iδ = IntervalMatrix(Diagonal(fill(IA.interval(δ), n)))
    else
        Iδ = Matrix(δ * I, n, n)
    end

    IδW = Iδ + 1 / 2 * δ^2 * A + 1 / 6 * δ^3 * A²
    M = IδW

    if order > 2
        # i = 2
        αᵢ₊₁ = 6 # factorial of (i+1)
        Aⁱ = A²
        δⁱ⁺¹ = δ^3
        @inbounds for i in 3:order
            αᵢ₊₁ *= i + 1
            δⁱ⁺¹ *= δ
            Aⁱ *= A
            M += (δⁱ⁺¹ / αᵢ₊₁) * Aⁱ
        end
    end
    E = _exp_remainder(A, δ, order; n=n)
    return M + E * δ
end

end  # module
