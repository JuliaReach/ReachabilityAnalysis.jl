using ..Overapproximate: _convert_or_overapproximate

"""
    GLGM06{N, AM, S, D, NG, P, RM} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
linear systems using zonotopes.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `FirstOrderZonotope`) approximation
                        model; see `Notes` below for possible options
- `max_order`        -- (optional, default: `5`) maximum zonotope order
- `static`           -- (optional, default: `false`) if `true`, convert the problem data
                        to statically sized arrays
- `dim`              -- (optional default: `missing`) ambient dimension
- `ngens`            -- (optional, default: `missing`) number of generators
- `preallocate`      -- (optional, default: `true`) if `true`, use the implementation
                        which preallocates the zonotopes prior to applying the update rule
- `reduction_method`    -- (optional, default: `GIR05()`) zonotope order reduction method used
- `disjointness_method` -- (optional, default: `NoEnclosure()`) method to check
                           disjointness between the reach-set and the invariant

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model
- `S`  -- value type associated to the `static` option
- `D`  -- value type associated to the dimension of the system
- `NG` -- value type associated to the number of generators
- `P`  -- value type associated to the `preallocate` option
- `RM` -- type associated to the reduction method

The sole parameter which doesn't have a default value is the step-size,
associated to the type parameter `N`. Parameters `D` and `NG` are optionally
specified (default to `Missing`). These parameters are needed for implementations
that require the size of the zonotopes to be known (fixed) at compile time, namely
the `static=true` version of this algorithm. Otherwise, the number of generators
is not necessarily fixed.

The default approximation model is

```julia
approx_model=FirstOrderZonotope()
```
Here, `FirstOrderZonotope` refers to the forward-time adaptation of the approximation model
from Lemma 3 in [[FRE11]](@ref). Some of the options to compute this approximation can be specified,
see the documentation of `FirstOrderZonotope` for details.

## References

The main ideas behind this algorithm can be found in [[GIR05]](@ref) and [[GLGM06]](@ref).
These methods are discussed at length in the dissertation [[LG09]](@ref).

Regarding the zonotope order reduction methods, we refer to [[COMB03]](@ref),
[[GIR05]](@ref) and the review article [[YS18]](@ref).

Regarding the approximation model, we use an adaptation of a result in [[FRE11]](@ref).
"""
struct GLGM06{N,AM,S,D,NG,P,RM,DM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
    static::S
    dim::D
    ngens::NG
    preallocate::P
    reduction_method::RM
    disjointness_method::DM
end

# TODO review setops "zonotope" / "lazy" options, used or ignored

# convenience constructor using symbols
function GLGM06(; δ::N,
                approx_model::AM=FirstOrderZonotope(),
                #CorrectionHull(order=10),
                max_order::Int=5,
                static::Bool=false,
                dim::Union{Int,Missing}=missing,
                ngens::Union{Int,Missing}=missing,
                preallocate::Bool=true,
                reduction_method::RM=LazySets.GIR05(),
                disjointness_method::DM=NoEnclosure()) where {N,AM,RM<:AbstractReductionMethod,
                                                              DM<:AbstractDisjointnessMethod}

    # algorithm with "preallocation" is only defined for the non-static case
    preallocate = !static
    n = ismissing(dim) ? missing : Val(dim)
    p = ismissing(ngens) ? missing : Val(ngens)
    return GLGM06(δ, approx_model, max_order, Val(static), n, p,
                  Val(preallocate), reduction_method, disjointness_method)
end

step_size(alg::GLGM06) = alg.δ
numtype(::GLGM06{N}) where {N} = N

@inline function setrep(alg::GLGM06{N,AM,Val{false}}) where {N,AM}
    return ST = Zonotope{N,Vector{N},Matrix{N}}
end

@inline function rsetrep(alg::GLGM06{N,AM,Val{false}}) where {N,AM}
    ST = setrep(alg)
    return RT = ReachSet{N,ST}
end

@inline function setrep(alg::GLGM06{N,AM,Val{true},Val{n},Val{p}}) where {N,AM,n,p}
    VT = SVector{n,N}
    MT = SMatrix{n,p,N,n * p}
    return ZT = Zonotope{N,VT,MT}
end

@inline function rsetrep(alg::GLGM06{N,AM,Val{true},Val{n},Val{p}}) where {N,AM,n,p}
    ZT = setrep(alg)
    return RT = ReachSet{N,ZT}
end

@inline function rsetrep(alg::GLGM06{N,AM,Val{true},Missing,Missing}) where {N,AM}
    throw(ArgumentError("the reach-set representation for a statically sized zonotope " *
                        "that the system's dimension argument has been defined; add `dim=...`" *
                        "to the algorithm constructor"))
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
#include("check.jl")
