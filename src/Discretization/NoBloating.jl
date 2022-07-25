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
    print(io, "`NoBloating` approximation model with:\n")
    print(io, "    - exponentiation method: $(alg.exp)\n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - invertibility assumption: $(alg.inv)\n")
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
