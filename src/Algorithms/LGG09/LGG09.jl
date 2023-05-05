using LazySets.Approximations: AbstractDirections

"""
    LGG09{N, AM, TN<:AbstractDirections} <: AbstractContinuousPost

Implementation of Girard - Le Guernic algorithm for reachability analysis
using support functions.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `Forward`) approximation model;
                    see `Notes` below for possible options
- `template`     -- (alias: `dirs`) struct that holds the directions (either lazily or concretely)
                    for each support function evaluation defining the template
- `vars`         -- (optional, default: all variables are computed) an integer or a vector of integers
                    specifying the variables of interest to automatically construct a template
                    using canonical directions; requires that `n` (or `dim`) is specified as well
- `static`       -- (optional, default: `false`) if `true`, use statically sized arrays
- `threaded`     -- (optional, default: `false`) if `true`, use multi-threading
                    to compute different template directions in parallel
- `sparse`       -- (optional, default: `false`) if `true`, convert the matrix exponential
                    obtained after discretization to a sparse matrix
- `cache`        -- (optional, default: `true`) if `true`, use a cache for intermediate
                    computations in the set recurrence loop
- `vars`         -- (optional, default: `missing`) used to specify the variables instead
                    of passing the template

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
struct LGG09{N,AM,VN,TN<:AbstractDirections{N,VN},S,VA} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    template::TN
    static::S
    threaded::Bool
    sparse::Bool
    cache::Bool
    vars::VA
end

# convenience constructor using symbols
function LGG09(; δ::N,
               approx_model::AM=Forward(; sih=:concrete, exp=:base, setops=:lazy),
               template=nothing,
               dirs=nothing, # alias for template
               vars=missing, # shortcut to specify variables of interest
               n=missing,    # required to construct the template from vars
               dim=missing,  # alias for n
               static::Bool=false,
               threaded::Bool=false,
               sparse::Bool=false,
               cache::Bool=true) where {N,AM}
    got_vars = !ismissing(vars)

    _n = missing
    if !ismissing(n)
        _n = n
    elseif !ismissing(dim)
        _n = dim
    end
    got_dim = !ismissing(_n)

    _template = missing
    if !isnothing(dirs)
        _template = dirs
    elseif !isnothing(template)
        _template = template
    end
    got_template = !ismissing(_template)

    if got_vars
        if !got_dim
            throw(ArgumentError("the ambient dimension, `n` (or `dim`), should be specified"))
        end
        directions = _get_template(vars, _n)
    elseif got_template
        if got_dim
            directions = _get_template(_template, _n)
        else
            directions = _get_template(_template)
        end
    else
        throw(ArgumentError("either `vars`, `dirs` or `template` should be specified"))
    end

    return LGG09(δ, approx_model, directions, Val(static), threaded, sparse, cache, vars)
end

_get_template(template::AbstractDirections) = template
_get_template(template::AbstractVector{N}) where {N<:Number} = CustomDirections([template])
function _get_template(template::AbstractVector{VT}) where {N<:Number,VT<:AbstractVector{N}}
    return CustomDirections(template)
end

_get_template(template::Symbol, n::Int) = _get_template(Val(template), n)
_get_template(::Val{:box}, n) = BoxDirections(n)
_get_template(::Val{:oct}, n) = OctDirections(n)

_get_template(vars::Int, n::Int) = _get_template((vars,), n)

function _get_template(vars::VecOrTuple, n::Int)
    m = length(vars)
    @assert m ≤ n "the number of variables should not exceed `n`"
    dirs = Vector{SingleEntryVector{Float64}}(undef, 2 * m)

    @inbounds for (i, vi) in enumerate(vars)
        d₊ = SingleEntryVector(vi, n, 1.0)
        d₋ = SingleEntryVector(vi, n, -1.0)
        dirs[i] = d₊
        dirs[i + m] = d₋
    end
    return CustomDirections(dirs, n)
end

step_size(alg::LGG09) = alg.δ
numtype(::LGG09{N}) where {N} = N
setrep(alg::LGG09) = setrep(alg.template)

function rsetrep(alg::LGG09{N,AM,VN,TN}) where {N,AM,VN,TN}
    SN = SubArray{N,1,Matrix{N},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
    return TemplateReachSet{N,VN,TN,SN}
end

include("reach_homog.jl")
include("reach_inhomog.jl")
include("post.jl")
