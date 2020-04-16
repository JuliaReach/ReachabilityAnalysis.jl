using LazySets.Approximations: AbstractDirections

"""
    LGG09{N, AM} <: AbstractContinuousPost

Implementation of Girard - Le Guernic algorithm for reachability analysis
using support functions.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional) approximation model for the discretization of the
                    ODE; see `Notes` below for details

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model

## References

See [1] and [2].

[1] Girard, A. (2005, March). Reachability of uncertain linear systems using zonotopes.
    In International Workshop on Hybrid Systems: Computation and Control (pp. 291-305).
    Springer, Berlin, Heidelberg.
    [Link](http://www-ljk.imag.fr/membres/Antoine.Girard/Publications/hscc2005.pdf)

[2] Girard, A., Le Guernic, C., & Maler, O. (2006, March).
    Efficient computation of reachable sets of linear time-invariant systems with inputs.
    In International Workshop on Hybrid Systems: Computation and Control (pp. 257-271).
    Springer, Berlin, Heidelberg.
    [Link](http://www-verimag.imag.fr/~maler/Papers/zonotope.pdf).
"""
@with_kw struct LGG09{N, AM, AD<:AbstractDirections} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
    template::AD
end

step_size(alg::LGG09) = alg.δ
numtype(::LGG09{N}) where {N} = N

include("reach.jl")
include("post.jl")
