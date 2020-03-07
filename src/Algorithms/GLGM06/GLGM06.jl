export GLGM06

"""
    GLGM06 <: AbstractContinuousPost

## Fields


- `δ`                   -- (optional, default: `1e-2`) step-size of the discretization
- `approximation_model` -- (optional, default `ForwardApproximation`) approximation model
                           for the discretization of the ODE
- `max_order`           -- (optional, default: `10`) maximum zonotope order

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
@with_kw struct GLGM06 <: AbstractContinuousPost
    δ::Float64
    approx_model::AbstractApproximationModel=ForwardApproximation(sih_method="concrete",
                                                                  exp_method="base",
                                                                  set_operations="zonotope",
                                                                  phi2_method="base")
    max_order::Int=10
    #setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
end

step_size(alg::GLGM06) = alg.δ
approx_model(alg::GLGM06) = alg.approx_model
max_order(alg::GLGM06) = alg.max_order

include("post.jl")
include("reach.jl")

function _project(R::ReachSet{N, ZT}, vars) where {N, ZT<:LazySets.AbstractZonotope}
    M = LazySets.Arrays.projection_matrix(vars, dim(R), N)
    return ReachSet(linear_map(M, set(R)), tspan(R))
end
