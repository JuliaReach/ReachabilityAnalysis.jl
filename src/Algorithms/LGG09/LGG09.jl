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
- `threaded`     -- (optional, default: `false`) if `true`, use multi-threading
                    to compute different template directions in parallel
- `sparse`       -- (optional, default: `false`) if `true`, convert the matrix exponential
                    obtained after discretization to a sparse matrix
- `cache`        -- (optional, default: `true`) if `true`, use a cache for intermediate
                    computations in the set recurrence loop

## Notes

The type fields are:

- `N`        -- number type of the step-size
- `AM`       -- type of the approximation model
- `TN`       -- type of the abstract directions that define the template

The flag `threaded=true` specifies that the support functions
along different directions should be computed in parallel.
See the section on [Multi-threading](@ref) for details on how to setup
the number of threads.

## References

This is an implementation of the algorithm from [[LGG09]](@ref).

These methods are described at length in the dissertation [[LG09]](@ref).
"""
struct LGG09{N, AM, VN, TN<:AbstractDirections{N, VN}, S} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    template::TN
    static::S
    threaded::Bool
    sparse::Bool
    cache::Bool
end

# convenience constructor using symbols
function LGG09(; δ::N,
               approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:lazy),
               template,
               static::Bool=false,
               threaded::Bool=false,
               sparse::Bool=false,
               cache::Bool=true) where {N, AM}

    return LGG09(δ, approx_model, _get_template(template), Val(static), threaded, sparse, cache)
end

_get_template(template::AbstractDirections) = template
_get_template(template::AbstractVector{N}) where {N<:Number} = CustomDirections([template])
_get_template(template::AbstractVector{VT}) where {N<:Number, VT<:AbstractVector{N}} = CustomDirections(template)

step_size(alg::LGG09) = alg.δ
numtype(::LGG09{N}) where {N} = N
setrep(alg::LGG09) = setrep(alg.template)

function rsetrep(alg::LGG09{N, AM, VN, TN}) where {N, AM, VN, TN}
    SN = SubArray{N, 1, Matrix{N}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}
    return TemplateReachSet{N, VN, TN, SN}
end

include("reach_homog.jl")
include("reach_inhomog.jl")
include("post.jl")
