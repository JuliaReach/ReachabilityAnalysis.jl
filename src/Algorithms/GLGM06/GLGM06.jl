export GLGM06

"""
    GLGM06{N, AM<:AbstractApproximationModel} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
uncertain linear systems using zonotopes.

## Fields

- `δ`           -- step-size of the discretization
- `appro_model` -- (optional, default `_DEFAULT_APPROX_MODEL_GLGM06`) approximation model
                   for the discretization of the ODE; see `Notes` below
- `max_order`   -- (optional, default: `10`) maximum zonotope order

## Notes


## References

TODO: move these references to the general references

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
@with_kw struct GLGM06{N, AM<:AbstractApproximationModel} <: AbstractContinuousPost
    δ::N
    approx_model::AM=_DEFAULT_APPROX_MODEL_GLGM06
    max_order::Int=10
end

const _DEFAULT_APPROX_MODEL_GLGM06 = ForwardApproximation(sih_method="concrete", exp_method="base",
                                                          set_operations="zonotope", phi2_method="base") # TODO use Val{...}

include("post.jl")
include("reach.jl")

#=
struct GLGM06{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
end

function DEFAULT_APPROX_MODEL

GLGM06(; δ, approx_model=DEFAULT_APPROX_MODEL(::GLGM06), max_order=10)
    = GLGM06{typeof(δ), ...)
=#
