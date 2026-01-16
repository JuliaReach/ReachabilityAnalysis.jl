"""
    A20{N} <: AbstractContinuousPost

Implementation of the reachability method for large linear systems with uncertain
inputs in the Krylov subspace from [Althoff20](@citet).

### Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `FirstOrderZonotope()`) approximation model
- `max_order`    -- (optional, default: `5`) maximum zonotope order

### References

See [Althoff20](@citet) and references therein.
"""
struct A20{N,AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
end

# convenience constructor using keywords
# TODO change `FirstOrderZonotope` default
function A20(; δ::N, approx_model::AM=FirstOrderZonotope(), max_order::Int=5) where {N,AM}
    return A20(δ, approx_model, max_order)
end

step_size(alg::A20) = alg.δ
numtype(::A20{N}) where {N} = N

function rsetrep(::A20{N}) where {N}
    VT = Vector{N}
    MT = Matrix{N}
    return ReachSet{N,Zonotope{N,VT,MT}}
end
