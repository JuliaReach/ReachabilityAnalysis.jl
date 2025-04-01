# ===============================
# First-order approximation model
# ===============================
module FirstOrderModule

using ..DiscretizationModule
using ..Exponentiation: _exp, BaseExp
using LinearAlgebra
using MathematicalSystems
using LazySets
using Reexport
using ReachabilityBase.Comparison: isapproxzero

export FirstOrder

@reexport import ..DiscretizationModule: discretize

const CLCS = ConstrainedLinearContinuousSystem
const CLCCS = ConstrainedLinearControlContinuousSystem

"""
    FirstOrder{EM} <: AbstractApproximationModel

First-order approximation model.

### Fields

- `exp` -- (optional, default: `BaseExp`) exponentiation method

### Algorithm

The transformations are [LeGuernicG10](@cite):

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← CH(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU ⊕ B_ε)``

where ``B_ε`` is the input ball of radius ``ε`` centered in the origin.
"""
struct FirstOrder{EM} <: AbstractApproximationModel
    exp::EM
end

# convenience constructor
function FirstOrder(; exp=BaseExp)
    return FirstOrder(exp)
end

function Base.show(io::IO, alg::FirstOrder)
    print(io, "`FirstOrder` approximation model with:\n")
    print(io, "    - exponentiation method: $(alg.exp)\n")
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", alg::FirstOrder) = print(io, alg)

# -----------------------------------------------
# FirstOrder approximation: Homogeneous case
# -----------------------------------------------

function discretize(ivp::IVP{<:CLCS,<:AbstractZonotope}, δ, alg::FirstOrder)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)
    norm_A = opnorm(A, Inf)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # compute Ω0
    Xδ = Φ * X0
    α = (exp(norm_A * δ) - 1 - norm_A * δ) * norm(X0, Inf)
    Ω0 = ConvexHull(X0, Xδ + BallInf(zeros(n), α))

    # create result
    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# -------------------------------------------------
# FirstOrder approximation: Inhomogeneous case
# -------------------------------------------------

function discretize(ivp::IVP{<:CLCCS,<:AbstractZonotope}, δ, alg::FirstOrder)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)
    norm_A = opnorm(A, Inf)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # compute Ω0
    Xδ = Φ * X0
    U = next_set(inputset(ivp), 1)
    factor = (exp(norm_A * δ) - 1 - norm_A * δ)
    norm_U_over_A = isapproxzero(norm_A) ? 0 : norm(U, Inf) / norm_A
    norm_U_over_A = norm(U, Inf) / norm_A
    α = factor * (norm(X0, Inf) + norm_U_over_A)
    Ω0 = ConvexHull(X0, Xδ + δ * U + BallInf(zeros(n), α))

    # discretize inputs
    β = factor * norm_U_over_A
    V = δ * U + BallInf(zeros(n), β)

    # create result
    B = IdentityMultiple(one(eltype(A)), size(A, 1))
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, B, X, V)
    return InitialValueProblem(Sdis, Ω0)
end

end  # module
