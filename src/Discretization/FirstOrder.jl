# ===============================
# First-order approximation model
# ===============================

using ..Exponentiation: _exp

"""
    FirstOrder{EM} <: AbstractApproximationModel

First-order approximation model.

### Fields

- `exp` -- (optional, default: `BaseExp`) exponentiation method

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← CH(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU ⊕ B_ε)``

where ``B_ε`` is the input ball of radius ``ε`` centered in the origin.

### Reference

C. Le Guernic and A. Girard: Reachability analysis of linear systems using
support functions. NAHS 2010.
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
    return print(io, "    - exponentiation method: $(alg.exp)\n")
end

Base.show(io::IO, m::MIME"text/plain", alg::FirstOrder) = print(io, alg)

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
