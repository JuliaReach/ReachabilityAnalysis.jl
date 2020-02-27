export GLGM06

"""
    GLGM06 <: AbstractContinuousPost

## Fields


- `δ`                   -- (optional, default: `1e-2`) step-size of the discretization
- `approximation_model` -- (optional, default `ForwardApproximation`) approximation model;
                            valid options are:

    - `sih_method`     -- (optional, default: `"concrete"`)
    - `exp_method`     -- (optional, default: `"base"`)
    - `exp_method`     -- (optional, default: `"base"`)

- `max_order`     -- (optional, default: `10`) maximum zonotope order

## References

TODO: move to general references

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
    δ::Float64
    approximation_model::AbstractApproximationModel=ForwardApproximation(sih_method="concrete",
                                exp_method="base", set_operations="zonotope", phi2_method="base")
    max_order::Int=10
end

step_size(alg::GLGM06) = alg.δ
approximation_model(alg::GLGM06) = alg.approximation_model
max_order(alg::GLGM06) = alg.max_order

include("post.jl")
include("reach.jl")
#include("project.jl")
