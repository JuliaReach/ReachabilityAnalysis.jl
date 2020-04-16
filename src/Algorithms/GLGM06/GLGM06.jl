"""
    GLGM06{N, AM} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
uncertain linear systems using zonotopes.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default `_DEFAULT_APPROX_MODEL_GLGM06`) approximation model
                    for the discretization of the ODE; see `Notes` below
- `max_order`    -- (optional, default: `10`) maximum zonotope order

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model

## References

See [xxx] and [yyy]
"""
@with_kw struct GLGM06{N, AM} <: AbstractContinuousPost
    δ::N
    # TODO set_operations="zonotope" used or ignored
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
    max_order::Int=10
    static::Bool=false
end

step_size(alg::GLGM06) = alg.δ
numtype(::GLGM06{N}) where {N} = N

function rsetrep(alg::GLGM06{N}) where {N}
    if !alg.static
        RT = ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}
    else
        error("not implemented")
        # n : other parameter of GLGM06 / or just a field alg.dim
        # p : order * dimension eg. alg.max_order * alg.dim
        # VT = SVector{n, N}
        # MT = SMatrix{n, p, N, n*p}
        # ZT = Zonotope{N, VT, MT}
        # RT = ReachSet{N, ZT}
    end
    return RT
end

include("post.jl")
include("reach.jl")
#include("check.jl")
