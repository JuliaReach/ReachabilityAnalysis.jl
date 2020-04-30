using LazySets.Approximations: AbstractDirections

"""
    LGG09{N, AM, TN<:AbstractDirections} <: AbstractContinuousPost

Implementation of Girard - Le Guernic algorithm for reachability analysis
using support functions.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `Forward`) approximation model;
                    see `Notes` below for possible options
- `template`     -- struct that holds the directions (either lazily or concretely)
                    for each support function evaluation defining the template
- `static`       -- (optional, default: `false`) if `true`, use statically sized arrays
- `threaded`     -- (optional, default: `true`) if `true`, use multi-threading parallelism
                    to compute the support function along each direction

## Notes

The type fields are:

- `N`        -- number type of the step-size
- `AM`       -- type of the approximation model
- `TN`       -- type of the abstract directions that define the template

## References

The is an implementation of the algorithm from [[LGG09]](@ref).

These methods are described at length in the dissertation [[LG09]](@ref).
"""
struct LGG09{N, AM, VN, TN<:AbstractDirections{N, VN}, S} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    template::TN
    static::S
    threaded::Bool
end

# convenience constructor using symbols
function LGG09(; δ::N,
               approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:lazy),
               template::TN,
               static::Bool=false,
               threaded::Bool=true) where {N, AM, TN}
    return LGG09(δ, approx_model, template, Val(static), threaded)
end

step_size(alg::LGG09) = alg.δ
numtype(::LGG09{N}) where {N} = N

function rsetrep(alg::LGG09{N, AM, TN}) where {N, AM, TN}
    RT = TemplateReachSet{N, eltype(TN), TN, Vector{N}}
    # TODO: SN is also Vector{N} or SVector ? if alg.static = true ?
    return RT
end

include("reach_homog.jl")
include("reach_inhomog.jl")
include("post.jl")
