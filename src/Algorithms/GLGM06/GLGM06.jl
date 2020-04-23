"""
    GLGM06{N, AM, D, NG, RM} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
linear systems using zonotopes.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `Forward`) approximation model;
                        see `Notes` below for possible options
- `max_order`        -- (optional, default: `5`) maximum zonotope order
- `static`           -- (optional, default: `false`) if `true`, convert the problem data
                        to statically sized arrays
- `dim`              -- (optional default: `missing`) ambient dimension
- `ngens`            -- (optional, default: `missing`) number of generators
- `preallocate`      -- (optional, default: `true`) if `true`, use the implementation
                        which preallocates the zonotopes prior to applying the update rule
- `reduction_method` -- (optional, default: `GIR05()`) zonotope order reduction method used

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model
- `D`  -- refers to the dimension of the system
- `NG` -- refers to the number of generators
- `RM` -- type associated to the reduction method

The sole parameter which doesn't have a default value is the step-size,
associated to the type parameter `N`. Parameters `D` and `NG` are optionally
specified (default to `Missing`). These parameters are needed for implementations
that require the size of the zonotopes to be known (fixed) at compile time, namely
the `static=true` version of this algorithm. Otherwise, the number of generators
is not necessarily fixed.

The default approximation model is

```julia
approx_model=Forward(sih=:concrete, exp=:base, setops=:lazy)
```
Here, `Forward` refers to the forward-time adaptation of the approximation model
from Lemma 3 in [[FRE11]](@ref). Some of the options to compute this approximation can be specified,
see the documentation of `Forward` for details.

## References

The main ideas behind this algorithm can be found in [[GIR05]](@ref) and [[GLGM06]](@ref).
These methods are discussed at length in the dissertation [[LG09]](@ref).

Regarding the zonotope order reduction methods, we refer to [[COMB03]](@ref),
[[GIR05]](@ref) and the review article [[YS18]](@ref).

Regarding the approximation model, we use an adaptation of a result in [[FRE11]](@ref).
"""
struct GLGM06{N, AM, S, D, NG, P, RM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:lazy)     # TODO review setops "zonotope" / "lazy" options, used or ignored
    max_order::Int
    static::S
    dim::D
    ngens::NG
    preallocate::P
    reduction_method::RM
end

# convenience constructor using symbols
function GLGM06(; δ::N,
               approx_model::AM=CorrectionHull(order=10, exp=:base),
               max_order::Int=5,
               static::Bool=false,
               dim::Union{Int, Missing}=missing,
               ngens::Union{Int, Missing}=missing,
               preallocate::Bool=true,
               reduction_method::RM=GIR05()) where {N, AM, RM}
    n = ismissing(dim) ? missing : Val(dim)
    p = ismissing(ngens) ? missing : Val(ngens)
    return GLGM06(δ, approx_model, max_order, Val(static),
                  Val(n), Val(p), Val(preallocate), reduction_method)
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
include("reach_homog.jl")
include("reach_inhomog.jl")
#include("check.jl")
