# ============================================================================
# Composition operation for one step forward in time with one step backward
# ============================================================================

"""
    StepIntersect{DM<:AbstractApproximationModel} <: AbstractApproximationModel

Approximation model that composes (intersecting) one step forward with one step backward.

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
struct StepIntersect{DM<:AbstractApproximationModel, SO} <: AbstractApproximationModel
    model::DM
    setops::SO
end

hasbackend(alg::StepIntersect) = !isnothing(alg.model)

# convenience constructor using symbols
function StepIntersect(; setops=:lazy)
    return StepIntersect(Forward(), _alias(setops))
end

$
function Base.show(io::IO, alg::Forward)
    print(io, "`StepIntersect` approximation model with: \n")
    print(io, "    - model: $(alg.model) \n")
    print(io, "    - set operations method: $(alg.setops)\n")
end

Base.show(io::IO, m::MIME"text/plain", alg::StepIntersect) = print(io, alg)

# ------------------------------------------------------------
# Homogeneous case
# ------------------------------------------------------------

# x' = Ax, x(0) ∈ X₀, x ∈ X 
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::StepIntersect)
    ivpd = discretize(ivp, δ, alg.model)
end
