"""
    ASB07{N, AM, RM} <: AbstractContinuousPost

Implementation of Althoff - Stursberg - Buss algorithm for reachability of
linear systems with uncertain parameters and inputs using zonotopes.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `Forward`) approximation model;
                        see `Notes` below for possible options
- `max_order`        -- (optional, default: `5`) maximum zonotope order
- `static`           -- (optional, default: `false`) if `true`, convert the problem data
                        to statically sized arrays
- `recursive`        -- (optional default: `true`) if `true`, use the implementation that
                        recursively computes each reach-set; otherwise, use the implementation
                        that unwraps the sequence until the initial set
- `reduction_method` -- (optional, default: `GIR05()`) zonotope order reduction method used

## Notes

The type fields are:

- `N`  -- number type of the step-size
- `AM` -- approximation model
- `RM` -- type associated to the reduction method

The sole parameter which doesn't have a default value is the step-size,
associated to the type parameter `N`.

The `static=true` version of this algorithm is not implemented yet.

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
@with_kw struct ASB07{N, AM, RM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=CorrectionHull(order=10, exp=:base)
    max_order::Int=5
    static::Bool=false
    recursive::Bool=true
    reduction_method::RM=GIR05()
end

step_size(alg::ASB07) = alg.δ
numtype(::ASB07{N}) where {N} = N
function rsetrep(alg::ASB07{N}) where {N}
    if !alg.static
        RT = ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}
    else
        error("not implemented yet")
        #=
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
        =#
    end
    return RT
end

include("post.jl")
include("reach_homog.jl")
