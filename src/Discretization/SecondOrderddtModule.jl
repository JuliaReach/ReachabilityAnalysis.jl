# ====================================
# Second-order approximation from d/dt
# ====================================
module SecondOrderddtModule

using ..DiscretizationModule
using ..Exponentiation: _exp, _alias, BaseExp, Φ₂
using ..ApplySetops: _apply_setops
using LinearAlgebra
using MathematicalSystems
using LazySets
using Reexport

export SecondOrderddt

@reexport import ..DiscretizationModule: discretize

const CLCS = ConstrainedLinearContinuousSystem

"""
    SecondOrderddt{EM, SO, SI, IT, BT} <: AbstractApproximationModel

Second-order approximation model used in the tool `d/dt`.
It can be used for overapproximation and underapproximation.

### Fields

- `oa`      -- (optional, default: `true`) flag to choose between
               overapproximation and underapproximation
- `exp`     -- (optional, default: `BaseExp`) exponentiation method
- `setops`  -- (optional, default: `:lazy`) set operations method
- `sih`     -- (optional, default: `:concrete`) way to compute the symmetric
               interval hull
- `inv`     -- (optional, default: `false`) if `true`, assume that the state
               matrix is invertible and use its inverse in the `Φ` functions
- `backend` -- (optional, default: `nothing`) used if the algorithm needs to
               apply concrete polyhedral computations

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← bloat(CH(\\mathcal{X}_0, Φ\\mathcal{X}_0), ε)``

where ``bloat(\\mathcal{X}, ε)`` bloats the set ``\\mathcal{X}`` with the value
``ε``. If `oa == false`, the bloating acts in an inverted way and shrinks the
set.

### Reference

E. Asarin, T. Dang, O. Maler, O. Bournez: *Approximate reachability analysis of
piecewise-linear dynamical systems*. HSCC 2000.
"""
struct SecondOrderddt{EM,SO,SI,IT,BT} <: AbstractApproximationModel
    oa::Bool
    exp::EM
    setops::SO
    sih::SI
    inv::IT
    backend::BT
end

hasbackend(alg::SecondOrderddt) = !isnothing(alg.backend)

# convenience constructor using symbols
function SecondOrderddt(; oa::Bool=true, exp=BaseExp, setops=:lazy, sih=:concrete,
                        inv=false, backend=nothing)
    return SecondOrderddt(oa, _alias(exp), _alias(setops), Val(sih), Val(inv), backend)
end

function Base.show(io::IO, alg::SecondOrderddt)
    print(io, "`SecondOrderddt` approximation model with:\n")
    print(io, "    - $(alg.oa ? "over" : "under")approximation\n")
    print(io, "    - exponentiation method: $(alg.exp)\n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - symmetric interval hull method: $(alg.sih)\n")
    print(io, "    - invertibility assumption: $(alg.inv)\n")
    print(io, "    - polyhedral computations backend: $(alg.backend)\n")
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", alg::SecondOrderddt) = print(io, alg)

# ------------------------------------------------------------
# SecondOrderddt Approximation: Homogeneous case
# ------------------------------------------------------------

# if A == |A|, then Φ can be reused in the computation of Φ₂(|A|, δ)
function discretize(ivp::IVP{<:CLCS,<:LazySet}, δ, alg::SecondOrderddt)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    Y = ConvexHull(X0, Φ * X0)
    Y = _apply_setops(Y, alg)

    # bloating
    ε = _estimate_bloating_value_ddt(A, X0, δ)
    if !alg.oa
        ε = -ε
    end
    Ω0 = Bloating(Y, ε)

    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

function _estimate_bloating_value_ddt(A, X0, δ)
    norm_A_δ = opnorm(A) * δ
    u = 1 / 8 * norm_A_δ^2
    v = exp(norm_A_δ) - 1 - norm_A_δ - (norm_A_δ^2 / 2)
    return (u + v) * norm(X0)
end

# outer bloating a polytope into a polytope
function _bloat(P::LazySet, ε)
    A, b = tosimplehrep(P)
    bε = copy(b)
    for i in eachindex(b)
        bε[i] += ε * opnorm(transpose(A[i, :]))
    end
    return HPolytope(A, bε)
end

end  # module
