# ==================================
# Forward approximation
# ==================================
module ForwardModule

using ..DiscretizationModule
using ..Exponentiation: _exp, _alias, BaseExp, elementwise_abs, Φ₂
using ..ApplySetops: _apply_setops
using IntervalArithmetic: hull
using MathematicalSystems: ConstrainedLinearContinuousSystem,
                           ConstrainedLinearControlContinuousSystem,
                           ConstrainedLinearControlDiscreteSystem,
                           ConstrainedLinearDiscreteSystem, IVP,
                           IdentityMultiple, initial_state, inputset,
                           state_matrix, stateset
using LazySets: ConvexHull, Interval, LazySet, symmetric_interval_hull, ⊕
using Reexport: @reexport

export Forward

@reexport import ..DiscretizationModule: discretize

const CLCS = ConstrainedLinearContinuousSystem
const CLCCS = ConstrainedLinearControlContinuousSystem

"""
    Forward{EM, SO, SI, IT, BT} <: AbstractApproximationModel

Forward approximation model.

### Fields

- `exp`     -- exponentiation method
- `setops`  -- set operations method
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

For the definition of the sets ``E_ψ`` and ``E^+`` see [FrehseGDCRLRGDM11](@cite).
The `Backward` method uses ``E^-``.
"""
struct Forward{EM,SO,SI,IT,BT} <: AbstractApproximationModel
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
    print(io, "    - invertibility assumption: $(alg.inv)\n")
    print(io, "    - polyhedral computations backend: $(alg.backend)\n")
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", alg::Forward) = print(io, alg)

# ------------------------------------------------------------
# Forward Approximation: Homogeneous case
# ------------------------------------------------------------

# if A == |A|, then Φ can be reused in the computation of Φ₂(|A|, δ)
function discretize(ivp::IVP{<:CLCS,<:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)
    E₊ = sih(P2A_abs * sih((A * A) * X0, alg.sih), alg.sih)

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ E₊)
    Ω0 = _apply_setops(Ω0, alg)

    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return IVP(Sdis, Ω0)
end

function discretize(ivp::IVP{<:CLCS,Interval{N}}, δ, ::Forward) where {N}
    A = state_matrix(ivp)
    @assert size(A, 1) == 1
    X0 = initial_state(ivp)

    a = A[1, 1]
    aδ = a * δ
    Φ = exp(aδ)

    # use inverse method
    @assert !iszero(a) "the given matrix should be invertible"

    # A_abs = abs(a)
    # a_sqr = a * a
    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    #E⁺ = (P2A_abs * a_sqr) * RA._symmetric_interval_hull(X0).dat

    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    E⁺ = (Φ - one(N) - aδ) * convert(Interval, symmetric_interval_hull(X0)).dat

    Ω0 = Interval(hull(X0.dat, Φ * X0.dat + E⁺))

    X = stateset(ivp)
    # the system constructor creates a matrix
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return IVP(Sdis, Ω0)
end

# ------------------------------------------------------------
# Forward Approximation: Inhomogeneous case
# ------------------------------------------------------------

# TODO : specialize, add option to compute the concrete linear map
function discretize(ivp::IVP{<:CLCCS,<:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = elementwise_abs(A)
    Φcache = sum(A) == abs(sum(A)) ? Φ : nothing
    P2A_abs = Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)

    # TODO outsource to Exponentiation module and merge with _Eplus
    E⁺ = sih(P2A_abs * sih((A * A) * X0, alg.sih), alg.sih)

    U = next_set(inputset(ivp), 1)
    Eψ0 = sih(P2A_abs * sih(A * U, alg.sih), alg.sih)

    # discretize inputs
    V = δ * U ⊕ Eψ0
    V = _apply_setops(V, alg.setops)

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ V ⊕ E⁺)
    Ω0 = _apply_setops(Ω0, alg.setops)

    # create result
    B = IdentityMultiple(one(eltype(A)), size(A, 1))
    X = stateset(ivp)
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, B, X, V)
    return IVP(Sdis, Ω0)
end

end  # module
