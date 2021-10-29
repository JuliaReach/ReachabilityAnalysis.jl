# ==================================
# Approximation model with zonotopes
# ==================================

"""
    FirstOrderZonotope{EM} <: AbstractApproximationModel

First order approximation model that works with zonotopes.

### Fields

- `exp` -- (optional, default: `BaseExp`) exponentiation method

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← bloat(zono(CH(\\mathcal{X}_0, Φ\\mathcal{X}_0)), α + β)``

where ``bloat(\\mathcal{X}, ε)`` bloats the set ``\\mathcal{X}`` with the value
``ε``, ``zono(·)`` overapproximates its argument with a zonotope, and ``α`` and
``β`` are factors computed for the homogeneous system and the inputs,
respectively.

### Reference

A. Girard: Reachability of uncertain linear systems using zonotopes. HSCC 2005.
"""
struct FirstOrderZonotope{EM} <: AbstractApproximationModel
    exp::EM
end

# convenience constructor
function FirstOrderZonotope(; exp=BaseExp)
    return FirstOrderZonotope(exp)
end

function Base.show(io::IO, alg::FirstOrderZonotope)
    print(io, "`FirstOrderZonotope` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
end

Base.show(io::IO, m::MIME"text/plain", alg::FirstOrderZonotope) = print(io, alg)

# -----------------------------------------------
# FirstOrderZonotope approximation: Homogeneous case
# -----------------------------------------------

function discretize(ivp::IVP{<:CLCS, <:AbstractZonotope}, δ, alg::FirstOrderZonotope)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)
    norm_A = opnorm(A, Inf)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # shared code
    Z, α = _discretize_zonotope(Φ, X0, alg, δ, norm_A)

    # compute bloating of Z
    Ω0 = minkowski_sum(Z, BallInf(zeros(n), α))

    # create result
    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# -------------------------------------------------
# FirstOrderZonotope approximation: Inhomogeneous case
# -------------------------------------------------

function discretize(ivp::IVP{<:CLCCS, <:AbstractZonotope}, δ, alg::FirstOrderZonotope)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)
    norm_A = opnorm(A, Inf)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # shared code
    Z, α = _discretize_zonotope(Φ, X0, alg, δ, norm_A)

    # compute bloating factor β
    U = next_set(inputset(ivp), 1)
    μ = norm(U, Inf)
    β = ((exp(norm_A * δ) - 1) * μ) / norm_A

    # compute bloating of Z
    Ω0 = minkowski_sum(Z, BallInf(zeros(n), α + β))

    # discretize inputs
    V = BallInf(zeros(n), β)

    # create result
    B = IdentityMultiple(one(eltype(A)), size(A, 1))
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, B, X, V)
    return InitialValueProblem(Sdis, Ω0)
end

# -----------
# Common code
# -----------

function _discretize_zonotope(Φ, X0, alg::FirstOrderZonotope, δ, norm_A)
    c = center(X0)
    G = genmat(X0)
    n, p = size(G)

    # compute X(δ)
    Xδ = linear_map(Φ, X0)

    # compute zonotope covering convex hull CH(X(0) ∪ X(δ))
    Z = overapproximate(CH(X0, Xδ), Zonotope)

    # compute bloating factor α
    α = (exp(norm_A * δ) - 1 - norm_A * δ) * norm(X0, Inf)

    return Z, α
end
