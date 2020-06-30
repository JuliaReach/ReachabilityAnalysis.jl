"""
    ASB07{N, AM, RM, S, R} <: AbstractContinuousPost

Implementation of Althoff - Stursberg - Buss algorithm for reachability of
linear systems with uncertain parameters and inputs using zonotopes.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `Forward`) approximation model;
                        see `Notes` below for possible options
- `max_order`        -- (optional, default: `5`) maximum zonotope order
- `reduction_method` -- (optional, default: `GIR05()`) zonotope order reduction method used
- `static`           -- (optional, default: `false`) if `true`, convert the problem data
                        to statically sized arrays
- `recursive`        -- (optional default: `true`) if `true`, use the implementation that
                        recursively computes each reach-set; otherwise, use the implementation
                        that unwraps the sequence until the initial set

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- type of the approximation model
- `RM` -- type associated to the reduction method
- `S`  -- value type associated to the `static` option
- `R`  -- value type associated to the `recursive` option

The sole parameter which doesn't have a default value is the step-size,
associated to the type parameter `N`.

The default approximation model is

```julia
approx_model=CorrectionHull(order=10, exp=:base)
```
Here, `CorrectionHull` refers to an implementation of the interval matrix
approximation method described in [[ASB07]](@ref). For technicalities on
interval matrix operations, we refer to the package `IntervalMatrices.jl`.

## References

The main ideas behind this algorithm can be found in [[ASB07]](@ref).
These methods are discussed at length in the dissertation [[ALT10]](@ref).

Regarding the zonotope order reduction methods, we refer to [[COMB03]](@ref),
[[GIR05]](@ref) and the review article [[YS18]](@ref).
"""
struct ASB07{N, AM, RM, S, R, D, NG} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
    reduction_method::RM
    static::S
    recursive::R
    dim::D
    ngens::NG
end

# convenience constructor using symbols
function ASB07(; δ::N,
               approx_model::AM=CorrectionHull(order=10, exp=:base),
               max_order::Int=5,
               reduction_method::RM=GIR05(),
               static::Bool=false,
               recursive::Bool=true,
               dim::Union{Int, Missing}=missing,
               ngens::Union{Int, Missing}=missing) where {N, AM, RM}

    n = ismissing(dim) ? missing : Val(dim)
    p = ismissing(ngens) ? missing : Val(ngens)
    return ASB07(δ, approx_model, max_order, reduction_method, Val(static), Val(recursive), n, p)
end

step_size(alg::ASB07) = alg.δ
numtype(::ASB07{N}) where {N} = N

function rsetrep(alg::ASB07{N, AM, RM, Val{false}}) where {N, AM, RM}
    RT = ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}
end

function rsetrep(alg::ASB07{N, AM, RM, S, R, Val{n}, Val{p}}) where {N, AM, RM, S, R, n, p}
    VT = SVector{n, N}
    MT = SMatrix{n, p, N, n*p}
    ZT = Zonotope{N, VT, MT}
    RT = ReachSet{N, ZT}
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
