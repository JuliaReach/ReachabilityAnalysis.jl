export GLGM06

"""
    GLGM06 <: AbstractContinuousPost

## Fields

- `δ`              -- (optional, default: `1e-2`) step-size of the discretization
- `discretization` -- (optional, default: `"forward"`)
- `sih_method`     -- (optional, default: `"concrete"`)
- `exp_method`     -- (optional, default: `"base"`)
- `max_order`      -- (optional, default: `10`)

## References

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
@with_kw struct GLGM06 <: AbstractContinuousPost
    δ::Float64=1e-2
    discretization::String="forward"
    sih_method::String="concrete"
    exp_method::String="base"
    max_order::Int=10
end

function init()
    # normalize
    # discretize
    # return system ready for calling the post operation post
end

#include("post.jl")
include("reach.jl")
#include("project.jl")
