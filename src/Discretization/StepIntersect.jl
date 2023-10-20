# ============================================================================
# Composition operation for one step of an approximation model with one step
# backward of the same model
# ============================================================================

using ..ExponentiationModule: _exp, _alias

"""
    StepIntersect{DM<:AbstractApproximationModel} <: AbstractApproximationModel

Approximation model that composes (intersecting) one step Forward of a given
model with one step backward of the same model.

### Fields

- `model`  -- approximation model
- `setops` -- set operations method

### Notes

Let ``x' = Ax`` with ``x(0) ∈ X₀``. This methods consists of:

- Compute the discretized system with step-size ``δ`` obtaining ``Ω0`` and the given
  approxatmion model `method`.

- Compute the (lazy) linear map ``ΦX₀``. This set constains of the (exact) reachable
  states at the time point ``δ``.

- Apply the approximation model `method` with initial condition ``ΦX₀`` one step
  backward in time, with the state transition matrix ``-A``, call this reach-set
  ``Ω₀₊``

- Intersect ``Ω₀₋`` and ``Ω₀₊`` and return such set. The intersection is done either
  lazily or concretely depending on the specified `setops` field.
"""
struct StepIntersect{DM<:AbstractApproximationModel,SO} <: AbstractApproximationModel
    model::DM
    setops::SO
end

hasbackend(alg::StepIntersect) = !isnothing(alg.model)

# convenience constructor using symbols
function StepIntersect(alg=Forward(); setops=alg.setops)
    return StepIntersect(alg, _alias(setops))
end

function Base.show(io::IO, alg::StepIntersect)
    print(io, "`StepIntersect` approximation model with:\n")
    print(io, "    - model: $(alg.model)\n")
    return print(io, "    - set operations method: $(alg.setops)\n")
end

Base.show(io::IO, m::MIME"text/plain", alg::StepIntersect) = print(io, alg)

# ------------------------------------------------------------
# Homogeneous case
# ------------------------------------------------------------

# x' = Ax, x(0) ∈ X₀, x ∈ X
function discretize(ivp::IVP{<:CLCS,<:LazySet}, δ, alg::StepIntersect)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)

    ivpd = discretize(ivp, δ, alg.model)
    Φ = state_matrix(ivpd)
    Ω0₊ = initial_state(ivpd)

    Aneg = -A
    ΦX0 = Φ * X0
    ivpneg = IVP(CLCS(Aneg, X), ΦX0)
    ivpnegd = discretize(ivpneg, δ, alg.model)
    Ω0₋ = initial_state(ivpnegd)

    Ω0 = _apply_setops(Ω0₊ ∩ Ω0₋, alg.setops)
    return IVP(CLDS(Φ, X), Ω0)
end

# ------------------------------------------------------------
# Inhomogeneous case
# ------------------------------------------------------------
