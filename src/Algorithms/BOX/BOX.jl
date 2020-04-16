"""
    BOX{N, AM} <: AbstractContinuousPost

Implementation of an integrator for linear systems usng box approximations.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional) approximation model for the discretization of the
                    continuous ODE; see the `Notes` below for available options

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model

The default approximation model used in this algorithm is:

```julia
Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
```
This algorithm overapproximates the discrete sequence ``X_{k+1} = ΦX_k ⊕ V_k``
by computing a tight hyperrectangular over-approximation of ``X_{k+1}`` at each
step ``k ∈ \\mathbb{N}``.
"""
@with_kw struct BOX{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
end

step_size(alg::BOX) = alg.δ
numtype(::BOX{N}) where {N} = N

include("post.jl")
include("reach.jl")
