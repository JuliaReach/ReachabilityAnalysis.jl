"""
    INT{N, AM} <: AbstractContinuousPost

Implementation of an interval-based integrator for one-dimensional systems.

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
Forward(sih=:concrete, exp=:base, phi2=:base, setops=:Interval)
```
In particular, the `setops=:Interval` flag specifies that intermediate computations
in the discretization are done using interval arithmetic.
"""
@with_kw struct INT{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:Interval)
end

step_size(alg::INT) = alg.δ

include("post.jl")
include("reach.jl")
