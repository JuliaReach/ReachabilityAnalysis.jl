# ==================================
# Backward approximation
# ==================================

using ..Exponentiation: _alias

"""
    Backward{EM, SO, SI, IT, BT} <: AbstractApproximationModel

Backward approximation model.

### Fields

- `exp`     -- exponentiation method
- `setops`  -- set opertaions method
- `sih`     -- symmetric interval hull
- `inv`     -- (optional, default: `false`) if `true`, assume that the state matrix
               is invertible and use its inverse in the `Φ` functions
- `backend` -- (optional, default: `nothing`) used if the algorithm needs to apply
               concrete polyhedral computations

### Algorithm

The transformations are:

- ``Φ ← \\exp(Aδ)``,
- ``Ω_0 ← CH(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU(0) ⊕ E_ψ(U(0), δ) ⊕ E^-(\\mathcal{X}_0, δ))``,
- ``V(k) ← δU(k) ⊕ E_ψ(U(k), δ)``.

Here we allow ``U`` to be a sequence of time varying non-deterministic input sets.

For the definition of the sets ``E_ψ`` and ``E^-`` see [[FRE11]](@ref).
The `Forward` method uses ``E^+``.
"""
struct Backward{EM,SO,SI,IT,BT} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
    inv::IT
    backend::BT
end

# convenience constructor using symbols
function Backward(; exp=BaseExp, setops=:lazy, sih=:concrete, inv=false, backend=nothing)
    return Backward(_alias(exp), _alias(setops), Val(sih), Val(inv), backend)
end

function Base.show(io::IO, alg::Backward)
    print(io, "`Backward` approximation model with:\n")
    print(io, "    - exponentiation method: $(alg.exp)\n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - symmetric interval hull method: $(alg.sih)\n")
    print(io, "    - invertibility assumption: $(alg.inv)\n")
    return print(io, "    - polyhedral computations backend: $(alg.backend)\n")
end

Base.show(io::IO, m::MIME"text/plain", alg::Backward) = print(io, alg)

# TODO: add corresponding `discrete` methods <<<<<
