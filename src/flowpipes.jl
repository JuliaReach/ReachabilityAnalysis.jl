export Flowpipe,
       HybridFlowpipe

# ================================
# Abstract types
# ================================

"""
    AbstractFlowpipe

Abstract type representing a flowpipe.

### Notes

A flowpipe is the set union of an array of reach-sets.
"""
abstract type AbstractFlowpipe end

"""
    project(F::AbstractFlowpipe, vars::AbstractVector)

Projects a flowpipe onto the subspace spanned by the given variables.

### Input

- `F`    -- flowpipe
- `vars` -- vector of variables

### Notes

This function can be used to project a high-dimensional flowpipe onto a
lower-dimensional space. The projection is lazy, and consists of mapping each set
`X` in the flowpipe to `MX`, where `M` is the projection matrix associated with
the given variables `vars`.
"""
function project(F::AbstractFlowpipe, vars::AbstractVector)
    return map(x -> LazySets.project(set(x), vars, LinearMap), array(F))
end

# ================================
# Flowpipes
# ================================

"""
    Flowpipe{ST, RT<:AbstractReachSet{ST}, VT<:AbstractVector{RT}} <: AbstractFlowpipe

Type that wraps a flowpipe.

### Fields

- `Xk`  -- set
- `ext` -- field used by extensions
"""
struct Flowpipe{T, RT<:AbstractReachSet{T}} <: AbstractFlowpipe
    Xk::Vector{RT}
    ext::Dict{Symbol, Any}
end

# constructor from empty ext dictionary
function Flowpipe(Xk::Vector{RT}) where {N, RT<:AbstractReachSet{N}}
    return Flowpipe(Xk, Dict{Symbol, Any}())
end

# TODO: use struct of array /// array of struct?
@inline array(fp::Flowpipe) = fp.Xk
LazySets.dim(fp::Flowpipe{ST, RT}) where {ST, RT<:AbstractReachSet{ST}} = dim(first(fp.Xk))
Base.length(fp::Flowpipe) = (length(fp.Xk),)

Base.first(fp::Flowpipe) = fp[1]
Base.last(fp::Flowpipe) = fp[end]
Base.firstindex(fp::Flowpipe) = 1
Base.lastindex(fp::Flowpipe) = length(fp.Xk)

# support indexing with ranges or with vectors of integers
Base.getindex(fp::Flowpipe, i::Int) = fp.Xk[i]
Base.getindex(fp::Flowpipe, i::Number) = fp[convert(Int, i)]
Base.getindex(fp::Flowpipe, I) = [fp[i] for i in I]

function tspan(fp::Flowpipe)
    tinit = inf(Δt(first(fp)))
    tend = sup(Δt(last(fp)))
    return IA.Interval{Float64}(tinit, tend)
end

# get the set of the flowpipe with the given index
function Base.getindex(fp::Flowpipe, t::Float64)
    1 <= i <= length(fp) || throw(BoundsError(fp, i))
    return fp.Xk[i]
end

# evaluate a flowpipe at a given time point: gives a reach set
# here it would be useful to layout the times contiguously in a vector ?
function (fp::Flowpipe)(t::Float64)
    for (i, X) in enumerate(fp.Xk)
        if t ∈ Δt(X) # exit on the first occurrence
            return fp[i]
        end
    end

    throw(ArgumentError("time $t does not belong to the time span, " *
            "$(tspan(fp)), of the given flowpipe"))
end

# evaluate a flowpipe at a given time interval: gives possibly more than one reach set
# ie. first and last sets and those in between them
function (fp::Flowpipe)(dt::IA.Interval{Float64})
    # here we assume that indices are of the form 1 .. n
    firstidx = 0
    lastidx = 0
    α = inf(dt)
    β = sup(dt)
    for (i, X) in enumerate(fp.Xk)
        if α ∈ Δt(X)
            firstidx = i
        end
        if β ∈ Δt(X)
            lastidx = i
        end
    end
    if firstidx == 0 || lastidx == 0
        throw(ArgumentError("time interval $Δt is not contained in the time span, " *
                "$(tspan(fp)), of the given flowpipe"))
    end
    return fp[firstidx:lastidx]
end

# ================================
# Hybrid flowpipe
# ================================

"""
    HybridFlowpipe{FT, VF<:AbstractVector{FT}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes.

### Fields

- `X`  -- set
- `Δt` -- time interval

### Notes

A `ReachSet` is a struct representing (an approximation of) the reachable states
for a given time interval. The type of the approximation is `ST`.
"""
struct HybridFlowpipe{FT, VF<:AbstractVector{FT}} <: AbstractFlowpipe
    Fk::VF
    ext::Dict{Symbol, Any}
end

@inline array(fp::HybridFlowpipe) = fp.Xk

#=
#dim(fp::Flowpipe{ST, RT}) where {ST, RT<:AbstractReachSet{ST}} = dim(first(fp.Xk))
# Base.getindex
=#

#=
# define the projection lazily?
# project(fp::Flowpipe, args; kwargs) -> lazy flowpipe with time
function project(fp::Flowpipe, args...; kwargs...)

    for X in fp.Xk # sets(fp.Xk)
        project(X, args...; kwargs...)
    end
end
=#

#=
"""
    project(Rsets, options)

Projects a sequence of sets according to the settings defined in the options.

### Input

- `Rsets`   -- solution of a reachability problem
- `options` -- options structure

### Notes

A projection matrix can be given in the options structure, or passed as a
dictionary entry.
"""
function project(Rsets::Vector{<:AbstractReachSet}, options::AbstractOptions)
    return project_reach(Rsets, options[:plot_vars], options)
end

TODO:

-> project a flowpipe

# helper function to take the cartesian product with the time variable
#fp::Flowpipe{FT, ST} where {FT, ST}


function add_time()
    flowpipe_with_time = Vector{ReachSet{CartesianProduct{Float64, IT64, ST}}}()
    add_time!()
end

function add_time!(F::Flowpipe{ST}) where {ST}
    @inbounds for X in fp
        Δt = X.Δt
        push!(fp, ReachSet(Δt × set(X), Δt))
    end
    return flowpipe_with_time
end
=#

#=
TESTS

using LazySets, Revise, ReachabilityAnalysis
N = Float64
X = Hyperrectangle(N[0, 0], N[1, 1])
R = ReachSet(X, 0 .. 1)
set(R) ==
tstart(R) == 0.0
tend(R) == 1.0
tspan(R) == 1.0
dim(R) == 2
overapproximate(project(X, [1]), Interval) == Interval(-1, 1)
=#
