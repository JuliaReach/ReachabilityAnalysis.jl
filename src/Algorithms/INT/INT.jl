"""
    INT{N, AM} <: AbstractContinuousPost

Implementation of a reachability method for linear one-dimensional systems
interval arithmetic.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `Forward`) approximation model;
                    see `Notes` below for possible options

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model

The default approximation model used in this algorithm is:

```julia
Forward(sih=:concrete, exp=:base, setops=:Interval)
```

In particular, the `setops=:Interval` flag specifies that intermediate computations
in the discretization are done using interval arithmetic. This allows for some
optimizations.

## References

This algorithm is essentially a non-decomposed version of the method in [BogomolovFFVPS18](@citet),
using intervals as set representation. For a general introduction we refer
to the dissertation [LeGuernicG09](@cite).

Regarding the approximation model, by default we use an adaptation of the method
presented in [FrehseGDCRLRGDM11](@citet).

Interval arithmetic operations are performed using the `IntervalArithmetic.jl`
package. Hence, the results are guaranteed to comply to the IEE754 standard with
respect to the floating-point operations using intervals.
"""
struct INT{N,AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
end

# convenience constructor using symbols
function INT(; δ::N,
             approx_model::AM=Forward(; sih=:concrete, exp=:base)) where {N,AM}
    return INT(δ, approx_model)
end

step_size(alg::INT) = alg.δ
numtype(::INT{N}) where {N} = N
setrep(::INT{N}) where {N} = Interval{N}
rsetrep(::INT{N}) where {N} = ReachSet{N,Interval{N}}

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
