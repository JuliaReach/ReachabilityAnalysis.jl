"""
    GLGM06{N, AM} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
uncertain linear systems using zonotopes.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default `_DEFAULT_APPROX_MODEL_GLGM06`) approximation model
                    for the discretization of the ODE; see `Notes` below
- `max_order`    -- (optional, default: `10`) maximum zonotope order
- `static`       -- (optional, default: `false`) if `true`, convert the problem data
                    to statically sized arrays
- `preallocate`  -- (optiona, default: `true`) if `true`, use the

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model

## References

See [xxx] and [yyy]
"""
@with_kw struct GLGM06{N, AM, D, NG} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy) # TODO set_operations="zonotope" used or ignored
    max_order::Int=5
    static::Bool=false
    dim::D=missing
    ngens::NG=missing
    preallocate::Bool=true
end

step_size(alg::GLGM06) = alg.δ
numtype(::GLGM06{N}) where {N} = N

function rsetrep(alg::GLGM06{N, AM, D}) where {N, AM, D}
    if !alg.static
        RT = ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}
    else
        @assert !ismissing(alg.dim) "the `static` option requires that the dimension " *
        "field of this algorithm is given, but it is $(alg.dim)"

        @assert !ismissing(alg.ngens) "the `static` option requires that the number of " *
        "generators is known, but it is $(alg.ngens)"

        n = alg.dim # dimension
        p = alg.ngens # number of generators
        VT = SVector{n, N}
        MT = SMatrix{n, p, N, n*p}
        ZT = Zonotope{N, VT, MT}
        RT = ReachSet{N, ZT}
    end
    return RT
end

include("post.jl")
include("reach.jl")
#include("check.jl")
