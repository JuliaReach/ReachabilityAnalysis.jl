# ============================================================
# Approximation model by optimization with template directions
# ============================================================

using JuMP

"""
    TemplateOptimization{DT, EM, PT<:Real, OT} <: AbstractApproximationModel

Approximation model that constructs an optimal convex approximation given a set
of template directions.

### Fields

- `D`   -- set of template directions
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
struct TemplateOptimization{DT, EM, PT<:Real, OT} <: AbstractApproximationModel
    D::DT
    exp::EM
    isbounded::Bool
    p::PT
    optimizer::OT
end

# convenience constructor
function TemplateOptimization(D; exp=BaseExp, isbounded=isbounding(D), p=Inf, optimizer)
    return TemplateOptimization(D, exp, isbounded, p, optimizer)
end

function Base.show(io::IO, alg::TemplateOptimization)
    print(io, "`TemplateOptimization` approximation model with: \n")
    print(io, "    - template directions: $(alg.D) \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - is the result known to be bounded? $(alg.isbounded) \n")
    print(io, "    - p-norm used: $(alg.p) \n")
    print(io, "    - optimizer: $(alg.optimizer) \n")
end

Base.show(io::IO, m::MIME"text/plain", alg::TemplateOptimization) = print(io, alg)

# ----------------------------------------------------
# TemplateOptimization approximation: Homogeneous case
# ----------------------------------------------------

function discretize(ivp::IVP{<:CLCS, <:AbstractZonotope}, δ, alg::TemplateOptimization)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # reachable states at time t
    R(t) = linear_map(exp(A * t), X0)

    # compute template approximation
    Ω0 = _discretize_template_optimization(δ, R, alg.D, alg.isbounded, alg.optimizer)

    # create result
    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# ------------------------------------------------------
# TemplateOptimization approximation: Inhomogeneous case
# ------------------------------------------------------

function discretize(ivp::IVP{<:CLCCS, <:AbstractZonotope}, δ, alg::TemplateOptimization)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)
    B = input_matrix(ivp)
    U = next_set(inputset(ivp), 1)

    # matrix exponential
    Φ = _exp(A, δ, alg.exp)

    # reachable states at time t
    R(t) = minkowski_sum(linear_map(exp(A * t), X0), linear_map(t * B, U))

    # compute template approximation
    Ω0 = _discretize_template_optimization(δ, R, alg.D, alg.isbounded, alg.optimizer)

    # discretize inputs
    # TODO make the specific method to obtain V an option
    # (currently we use one specific version from another approach)
    norm_A = opnorm(A, alg.p)
    μ = norm(U, alg.p)
    β = ((exp(norm_A * δ) - 1) * μ) / norm_A
    V = Ballp(alg.p, zeros(n), β)

    # create result
    B = IdentityMultiple(one(eltype(A)), size(A, 1))
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, B, X, V)
    return InitialValueProblem(Sdis, Ω0)
end

# -----------
# Common code
# -----------

function _discretize_template_optimization(δ, R, D, isbounded, optimizer)
    P = isbounded ? HPolytope() : HPolyhedron()
    for d in D
        f(t) = ρ(d, R(t))

        model = Model(typeof(optimizer))
        set_silent(model)
        @variable(model, 0 <= λ <= 1)
        register(model, :f, 1, f; autodiff = true)
        @NLobjective(model, Max, f(λ))
        optimize!(model)
        b = objective_value(model)

        c = HalfSpace(d, b)
        addconstraint!(P, c)
    end
    return P
end
