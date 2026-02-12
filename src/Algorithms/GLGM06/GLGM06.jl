"""
    GLGM06{N, AM, S, D, NG, P, RM} <: AbstractContinuousPost

Implementation of the Girard–Le Guernic–Maler algorithm for reachability of
linear systems using zonotopes.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `FirstOrderZonotope()`) approximation
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

The type parameters are:

- `N`  -- number type of the step-size
- `AM` -- approximation model
- `S`  -- value type of the `static` option
- `D`  -- value type of the dimension of the system `dim`
- `NG` -- value type of the number of generators `ngens`
- `P`  -- value type of the `preallocate` option
- `RM` -- type of the reduction method

The only parameter that does not have a default value is the step size `δ`,
associated with the type parameter `N`. Parameters `dim` and `ngens` are optionally
specified (default to `missing`). These parameters are needed for the cases
that require the size of the zonotopes to be known (fixed) at compile time, namely
the `static=true` version of this algorithm.

The default approximation model is [`FirstOrderZonotope`](@ref).

## References

The main ideas behind this algorithm can be found in [Girard05](@citet) and [GirardGM06](@citet).
These methods are discussed at length in the dissertation [LeGuernicG09](@cite).

Regarding the zonotope order reduction methods, we refer to [Combastel03, Girard05](@citet)
and the review article [YangS18](@cite).
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

# convenience constructor using symbols
function GLGM06(; δ::N,
                approx_model::AM=FirstOrderZonotope(),
                max_order::Int=5,
                static::Bool=false,
                dim::Union{Int,Missing}=missing,
                ngens::Union{Int,Missing}=missing,
                preallocate::Bool=true,
                reduction_method::RM=GIR05(),
                disjointness_method::DM=NoEnclosure()) where {N,AM,RM<:AbstractReductionMethod,
                                                              DM<:AbstractDisjointnessMethod}

    # algorithm with "preallocation" is only defined for the non-static case
    preallocate = preallocate && !static
    n = ismissing(dim) ? missing : Val(dim)
    p = ismissing(ngens) ? missing : Val(ngens)
    return GLGM06(δ, approx_model, max_order, Val(static), n, p,
                  Val(preallocate), reduction_method, disjointness_method)
end

step_size(alg::GLGM06) = alg.δ
numtype(::GLGM06{N}) where {N} = N

function setrep(::GLGM06{N,AM,Val{false}}) where {N,AM}
    return Zonotope{N,Vector{N},Matrix{N}}
end

function setrep(::GLGM06{N,AM,Val{true},Val{n},Val{p}}) where {N,AM,n,p}
    VT = SVector{n,N}
    MT = SMatrix{n,p,N,n * p}
    return Zonotope{N,VT,MT}
end

# TODO update error message (also requires `ngens`) and cover cases where only one argument is `missing`
function setrep(::GLGM06{N,AM,Val{true},Missing,Missing}) where {N,AM}
    throw(ArgumentError("the set representation for a statically sized zonotope " *
                        "that the system's dimension argument has been defined; add `dim=...`" *
                        "to the algorithm constructor"))
end

function rsetrep(alg::GLGM06{N}) where {N}
    ST = setrep(alg)
    return ReachSet{N,ST}
end
