using ..Overapproximate: _convert_or_overapproximate

"""
    A20{N} <: AbstractContinuousPost

Implementation of the reachability method for large linear systems with uncertain
inputs in the Krylov subspace from [Althoff20](@citet).

### Fields

- `δ`          -- step-size of the discretization
- `max_order`  -- (optional, default: `5`) maximum zonotope order

### References

See [Althoff20](@citet) and references therein.
"""
struct A20{N} <: AbstractContinuousPost
    δ::N
    max_order::Int
end

# convenience constructor using keywords
function A20(; δ::N, max_order::Int=5) where {N}
    return A20(δ, max_order)
end

step_size(alg::A20) = alg.δ
numtype(::A20{N}) where {N} = N

function rsetrep(::A20{N}) where {N}
    VT = Vector{N}
    MT = Matrix{N}
    return ReachSet{N,Zonotope{N,VT,MT}}
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
