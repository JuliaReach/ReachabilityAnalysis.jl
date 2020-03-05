export Flowpipe,
       HybridFlowpipe,
       flowpipe

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
    # TODO: add lazy option
    return map(x -> LazySets.project(set(x), vars, LinearMap), array(F))
end

# TODO: Implement LazySets interface; the flowpipe behaves like the union set array
#LazySets.ρ(d::AbstractVector, fp::AbstractFlowpipe) = ρ(d, set(R))
#LazySets.σ(d::AbstractVector, fp::AbstractFlowpipe) = σ(d, set(R))
#LazySets.dim(R::AbstractFlowpipe) = dim(set(R))

# ================================
# Flowpipes
# ================================

"""
    Flowpipe{N, RT<:AbstractReachSet{N}} <: AbstractFlowpipe

Type that wraps a flowpipe.

### Fields

- `Xk`  -- set
- `ext` -- extension dictionary; field used by extensions

### Notes

The dimension of the flowpipe corresponds to the dimension of the underlying
reach-sets; in this type, it is is assumed that the dimension is the same for
the different reach-sets.
"""
struct Flowpipe{N, RT<:AbstractReachSet{N}} <: AbstractFlowpipe
    Xk::Vector{RT}
    ext::Dict{Symbol, Any}
end

# getter functions
array(fp::Flowpipe) = fp.Xk
flowpipe(fp::Flowpipe) = fp

# constructor from empty extension dictionary
function Flowpipe(Xk::Vector{RT}) where {N, RT<:AbstractReachSet{N}}
    return Flowpipe(Xk, Dict{Symbol, Any}())
end

# iteration interface
Base.iterate(fp::Flowpipe) = iterate(fp.Xk)
Base.iterate(fp::Flowpipe, state) = iterate(fp.Xk, state)
Base.length(fp::Flowpipe) = length(fp.Xk)
Base.first(fp::Flowpipe) = fp[1]
Base.last(fp::Flowpipe) = fp[end]
Base.firstindex(fp::Flowpipe) = 1
Base.lastindex(fp::Flowpipe) = length(fp.Xk)

# abstract reach set interface
set(fp::Flowpipe) = throw(ArgumentError("to retrieve the array of sets represented by this flowpipe, " *
    "use the `array(...)` function, or use the function `set(...)` at a specific index, i.e. " *
    "`set(F[ind])`, or simply `set(F, ind)`, to get the reach-set with index `ind` of the flowpipe `F`"))
set(fp::Flowpipe, ind::Integer) = set(fp.Xk[ind])
setrep(fp::Flowpipe{N, RT}) where {N, RT} = setrep(RT)
setrep(::Type{Flowpipe{N, RT}}) where {N, RT} = setrep(RT)
@inline tstart(fp::Flowpipe) = tstart(first(fp))
@inline tend(fp::Flowpipe) = tend(last(fp))
@inline tspan(fp::Flowpipe) = IA.Interval(tstart(fp), tend(fp))
LazySets.dim(fp::Flowpipe) = dim(first(fp))

# support indexing with ranges or with vectors of integers
Base.getindex(fp::Flowpipe, i::Int) = fp.Xk[i]
Base.getindex(fp::Flowpipe, i::Number) = fp[convert(Int, i)]
Base.getindex(fp::Flowpipe, I::AbstractVector) = [fp[i] for i in I]

# get the set of the flowpipe with the given index
function Base.getindex(fp::Flowpipe, t::Float64)
    1 <= i <= length(fp) || throw(BoundsError(fp, i))
    return fp.Xk[i]
end

# evaluate a flowpipe at a given time point: gives a reach set
# here it would be useful to layout the times contiguously in a vector
# (see again array of struct vs struct of array)
function (fp::Flowpipe)(t::Float64)
    for (i, X) in enumerate(fp.Xk)
        if t ∈ tspan(X) # exit on the first occurrence
            return fp[i]
        end
    end
    throw(ArgumentError("time $t does not belong to the time span, " *
            "$(tspan(fp)), of the given flowpipe"))
end

# conversion from other numeric type
function (fp::Flowpipe)(t::Number)
    fp(Float64(t))
end

# evaluate a flowpipe at a given time interval: gives possibly more than one reach set
# i.e. first and last sets and those in between them
function (fp::Flowpipe)(dt::IA.Interval{Float64})
    # here we assume that indices are one-based, ie. form 1 .. n
    firstidx = 0
    lastidx = 0
    α = inf(dt)
    β = sup(dt)
    for (i, X) in enumerate(fp.Xk)
        if α ∈ tspan(X)
            firstidx = i
        end
        if β ∈ tspan(X)
            lastidx = i
        end
    end
    if firstidx == 0 || lastidx == 0
        throw(ArgumentError("the time interval $dt is not contained in the time span, " *
                "$(tspan(fp)), of the given flowpipe"))
    end
    return fp[firstidx:lastidx]
end

# ================================
# Hybrid flowpipe
# ================================

#=
# TODO: consider struct of array vs. array of struct
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
    # TODO: use a VectorOfArray{ST, 2, Vector{VT}}
    Fk::VF
    ext::Dict{Symbol, Any}
end
=#

#array(fp::HybridFlowpipe) = fp.Xk

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
