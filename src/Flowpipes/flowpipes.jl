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
    basetype(T::Type{<:AbstractFlowpipe})

Return the base type of the given flowpipe type (i.e., without type parameters).

### Input

- `T` -- flowpipe type, used for dispatch

### Output

The base type of `T`.
"""
basetype(T::Type{<:AbstractFlowpipe}) = Base.typename(T).wrapper

# LazySets interface: fallback behaves like UnionSetArray

"""
    LazySets.ρ(d::AbstractVector, fp::AbstractFlowpipe)

### Input

- `d`  -- direction
- `fp` -- flowpipe

### Output

The support function of the flowpipe along the given direction `d`.

### Notes

In this fallback implementation, the flowpipe behaves like the union of the
reach-sets, i.e. the implementation is analogue to that of a `LazySet.UnionSetArray`.
"""
function LazySets.ρ(d::AbstractVector, fp::AbstractFlowpipe)
    return map(Ri -> ρ(d, set(Ri)), array(fp)) |> maximum
end

"""
    LazySets.σ(d::AbstractVector, fp::AbstractFlowpipe)

### Input

- `d`  -- direction
- `fp` -- flowpipe

### Output

The support vector of the flowpipe along the given direction `d`.

### Notes

In this fallback implementation, the flowpipe behaves like the union of the
reach-sets, i.e. the implementation is analogue to that of a `LazySet.UnionSetArray`.
"""
function LazySets.σ(d::AbstractVector, fp::AbstractFlowpipe)
    return _σ_vec(d, array(fp))
end

"""
    LazySets.dim(fp::AbstractFlowpipe)

### Input

- `fp` -- flowpipe

### Output

An integer representing the ambient dimension of the flowpipe.
"""
function LazySets.dim(fp::AbstractFlowpipe)
    length(fp) > 0 || throw(ArgumentError("the dimension is not defined because this flowpipe is empty"))
    return dim(first(fp)) # assumes that the first set is representative
end

# iteration interface
@inline Base.iterate(fp::AbstractFlowpipe) = iterate(array(fp))
@inline Base.iterate(fp::AbstractFlowpipe, state) = iterate(array(fp), state)
@inline Base.length(fp::AbstractFlowpipe) = length(array(fp))
#@inline Base.size(fp::AbstractFlowpipe) = (length(array(fp)),)
@inline Base.first(fp::AbstractFlowpipe) = getindex(fp, 1)
@inline Base.last(fp::AbstractFlowpipe) = getindex(fp, lastindex(fp))
@inline Base.firstindex(fp::AbstractFlowpipe) = 1
@inline Base.lastindex(fp::AbstractFlowpipe) = length(array(fp))
@inline Base.eachindex(fp::AbstractFlowpipe) = eachindex(array(fp))

# support abstract reach set interface

set(fp::AbstractFlowpipe) = throw(ArgumentError("to retrieve the array of sets represented by this flowpipe, " *
    "use the `array(...)` function, or use the function `set(...)` at a specific index, i.e. " *
    "`set(F[ind])`, or simply `set(F, ind)`, to get the reach-set with index `ind` of the flowpipe `F`"))

"""
    set(fp::AbstractFlowpipe, ind::Integer)

Return the geometric set represented by this flowpipe at the given index.

## Input

- `fp`  -- flowpipe
- `ind` -- index (from `1` to `length(flowpipe)`)

## Output

The set wrapped by the flowpipe at the given index.
"""
set(fp::AbstractFlowpipe, ind::Integer) = set(getindex(array(fp), ind))

# time domain interface

"""
    tstart(fp::AbstractFlowpipe)

Return the initial time of this flowpipe.

### Input

- `fp` -- flowpipe

### Output

A float representing the initial time of the given flowpipe. The fallback is
computed by taking the initial time of the first reach-set.
"""
@inline tstart(fp::AbstractFlowpipe) = tstart(first(fp))

"""
    tend(fp::AbstractFlowpipe)

Return the final time of this flowpipe.

### Input

- `R` -- reach-set

### Output

A float representing the initial time of the given flowpipe. The fallback is
computed by taking the final time of the last reach-set.
"""
@inline tend(fp::AbstractFlowpipe) = tend(last(fp))

"""
    tspan(fp::AbstractFlowpipe)

Return time span of this flowpipe.

### Input

- `fp` -- flowpipe

### Output

The interval representing the time span of the given flowpipe. The fallback
is computed as `(tstart(fp), tend(fp))`, see `tstart(::AbstractFlowpipe)` and
`tend(::AbstractFlowpipe)` for details.
"""
@inline tspan(fp::AbstractFlowpipe) = TimeInterval(tstart(fp), tend(fp))

# assumes first set is representative
vars(fp::AbstractFlowpipe) = vars(first(fp))

# support indexing with ranges or with vectors of integers
# TODO add bounds checks?
Base.getindex(fp::AbstractFlowpipe, i::Int) = getindex(array(fp), i)
Base.getindex(fp::AbstractFlowpipe, i::Number) = getindex(array(fp), convert(Int, i))
Base.getindex(fp::AbstractFlowpipe, I::AbstractVector) = getindex(array(fp), I)

# get the set of the flowpipe with the given index
#function Base.getindex(fp::AbstractFlowpipe, t::Number)
    # annotate as a boundscheck
#    1 <= i <= length(fp) || throw(BoundsError(fp, i))
#    return getindex(fp, i)

#=
function Projection(fp::Flowpipe, vars::NTuple{D, T}) where {D, T<:Integer}

end
=#

#=
# inplace projection
function project!(fp::AbstractFlowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(fp)
    for X in Xk
        _project!(set(X), vars)
    end
    return fp
end
=#

# further setops
LazySets.is_intersection_empty(F::AbstractFlowpipe, Y::LazySet) where {N} = all(X -> _is_intersection_empty(X, Y), array(F))
Base.:⊆(F::AbstractFlowpipe, X::LazySet) = all(R ⊆ X for R in F)
Base.:⊆(F::AbstractFlowpipe, Y::AbstractLazyReachSet) = all(R ⊆ set(Y) for R in F)

# getter functions for hybrid systems
location(F::AbstractFlowpipe) = get(F.ext, :loc_id, missing)

# ================================
# Flowpipes
# ================================

"""
    Flowpipe{N, RT<:AbstractReachSet{N}, VRT<:AbstractVector{RT}} <: AbstractFlowpipe

Type that wraps a flowpipe.

### Fields

- `Xk`  -- set
- `ext` -- extension dictionary; field used by extensions

### Notes

The dimension of the flowpipe corresponds to the dimension of the underlying
reach-sets; in this type, it is is assumed that the dimension is the same for
the different reach-sets.
"""
struct Flowpipe{N, RT<:AbstractReachSet{N}, VRT<:AbstractVector{RT}} <: AbstractFlowpipe
    Xk::VRT
    ext::Dict{Symbol, Any}
end

# getter functions
@inline array(fp::Flowpipe) = fp.Xk
@inline flowpipe(fp::Flowpipe) = fp

# constructor from empty extension dictionary
function Flowpipe(Xk::AbstractVector{RT}) where {N, RT<:AbstractReachSet{N}}
    return Flowpipe(Xk, Dict{Symbol, Any}())
end

# undef initializer given a set type
function Flowpipe(::UndefInitializer, ST::Type{<:LazySet{N}}, k::Int) where {N}
    return Flowpipe(Vector{ReachSet{N, ST}}(undef, k), Dict{Symbol, Any}())
end

# undef initializer given a reach-set type
function Flowpipe(::UndefInitializer, RT::Type{<:AbstractReachSet}, k::Int)
    return Flowpipe(Vector{RT}(undef, k), Dict{Symbol, Any}())
end

# undef initializer given a continuous post-operator
function Flowpipe(::UndefInitializer, cpost::AbstractContinuousPost, k::Int)
    RT = rsetrep(cpost)
    return Flowpipe(Vector{RT}(undef, k), Dict{Symbol, Any}())
end

# constructor from a single reach-set
Flowpipe(R::AbstractReachSet) = Flowpipe([R])

function Base.similar(fp::Flowpipe{N, RT, VRT}) where {N, RT, VRT}
   return Flowpipe(VRT())
end

Base.IndexStyle(::Type{<:Flowpipe}) = IndexLinear()
Base.eltype(::Flowpipe{N, RT}) where {N, RT} = RT
Base.size(fp::Flowpipe) = (length(fp.Xk),)
Base.view(fp::Flowpipe, args...) = view(fp.Xk, args...)
Base.push!(F::Flowpipe, args...) = push!(F.Xk, args...)

numtype(::Flowpipe{N}) where {N} = N
setrep(fp::Flowpipe{N, RT}) where {N, RT} = setrep(RT)
setrep(::Type{<:Flowpipe{N, RT}}) where {N, RT} = setrep(RT)
rsetrep(fp::Flowpipe{N, RT}) where {N, RT} = RT
rsetrep(::Type{<:Flowpipe{N, RT}}) where {N, RT} = RT
numrsets(fp::Flowpipe) = length(fp)

function location(fp::Flowpipe)
    @assert haskey(fp.ext, :loc_id) "this flowpipe has not been assigned a location identifier"
    return fp.ext[:loc_id]
end

# evaluate a flowpipe at a given time point: gives a reach set
# here it would be useful to layout the times contiguously in a vector
# (see again array of struct vs struct of array)
function (fp::AbstractFlowpipe)(t::Number)
    Xk = array(fp)
    @inbounds for (i, X) in enumerate(Xk)
        if t ∈ tspan(X) # exit on the first occurrence
            if i < length(Xk) && t ∈ tspan(Xk[i+1])
                return view(Xk, i:i+1)
            else
                return X
            end
        end
    end
    throw(ArgumentError("time $t does not belong to the time span, " *
                        "$(tspan(fp)), of the given flowpipe"))
end

# evaluate a flowpipe at a given time interval: gives possibly more than one reach set
# i.e. first and last sets and those in between them
function (fp::Flowpipe)(dt::TimeInterval)
    # here we assume that indices are one-based, ie. form 1 .. n
    firstidx = 0
    lastidx = 0
    α = inf(dt)
    β = sup(dt)
    Xk = array(fp)
    for (i, X) in enumerate(Xk)
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
    return view(Xk, firstidx:lastidx)
end

function project(fp::Flowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(fp)
    # TODO: use projection of the reachsets
    if 0 ∈ vars # projection includes "time"
        # we shift the vars indices by one as we take the Cartesian prod with the time spans
        aux = vars .+ 1
        return map(X -> _project(convert(Interval, tspan(X)) × set(X), aux), Xk)
    else
        return map(X -> _project(set(X), vars), Xk) # TODO return Flowpipe ?
    end
end

project(fp::Flowpipe, vars::Int) = project(fp, (vars,))
project(fp::Flowpipe, vars::AbstractVector{<:Int}) = project(fp, Tuple(vars))
project(fp::Flowpipe; vars) = project(fp, Tuple(vars))
project(fp::Flowpipe, i::Int, vars) = project(fp[i], vars)

# concrete projection of a flowpipe for a given matrix
function project(fp::Flowpipe, M::AbstractMatrix; vars=nothing)
    Xk = array(fp)
    πfp = Flowpipe(map(X -> linear_map(M, X), Xk))
    if isnothing(vars)
        return πfp
    else
        return project(πfp, vars)
    end
end

# concrete projection of a flowpipe for a given direction
function project(fp::Flowpipe, dir::AbstractVector{<:AbstractFloat}; vars=nothing)
    Xk = array(fp)
    return Flowpipe(map(X -> project(X, dir, vars=vars), Xk))
end

# concrete linear map of a flowpipe for a given matrix
function linear_map(M::AbstractMatrix, fp::Flowpipe)
    Xk = array(fp)
    return Flowpipe(map(X -> linear_map(M, X), Xk))
end

"""
    shift(fp::Flowpipe{N, <:AbstractReachSet}, t0::Number) where {N}

Return the time-shifted flowpipe by the given number.

### Input

- `fp` -- flowpipe
- `t0` -- time shift

### Output

A new flowpipe such that the time-span of each constituent reach-set has been
shifted by `t0`.

### Notes

See also `Shift` for the lazy counterpart.
"""
function shift(fp::Flowpipe{N, <:AbstractReachSet}, t0::Number) where {N}
    return Flowpipe([shift(X, t0) for X in array(fp)], fp.ext)
end

"""
    convexify(fp::Flowpipe{N, <:AbstractLazyReachSet}) where {N}

Return a reach-set representing the convex hull array of the flowpipe.

### Input

- `fp` -- flowpipe

### Output

A reach-set that contains the convex hull array, `ConvexHullArray`, of the given
flowpipe.

### Notes

The time span of this reach-set is the same as the time-span of the flowpipe.

This function allocates an array to store the sets of the flowpipe.
"""
function convexify(fp::Flowpipe{N, <:AbstractLazyReachSet}) where {N}
    Y = ConvexHullArray([set(X) for X in array(fp)])
    return ReachSet(Y, tspan(fp))
end

"""
    convexify(fp::AbstractVector{<:AbstractLazyReachSet{N}}) where {N}

Return a reach-set representing the convex hull array of the array of the array of
reach-sets.

### Input

- `fp` -- array of reach-sets

### Output

A reach-set that contains the convex hull array, `ConvexHullArray`, of the given
flowpipe.

### Notes

The time span of this reach-set corresponds to the minimum (resp. maximum) of the
time span of each reach-set in `fp`.

This function allocates an array to store the sets of the flowpipe.
"""
function convexify(fp::AbstractVector{<:AbstractLazyReachSet{N}}) where {N}
    Y = ConvexHullArray([set(X) for X in fp])
    ti = minimum(tstart, fp)
    tf = minimum(tend, fp)
    return ReachSet(Y, TimeInterval(ti, tf))
end

# the dimension of sparse flowpipes is known in the type
function LazySets.dim(::Flowpipe{N, SparseReachSet{N, ST, D}}) where {N, ST, D}
    return D
end

tstart(F::Flowpipe, arr::UnitRange) = tstart(view(array(F), arr), contiguous=true)
tstart(F::Flowpipe, arr::AbstractVector) = tstart(view(array(F), arr))
tend(F::Flowpipe, arr::UnitRange) = tend(view(array(F), arr), contiguous=true)
tend(F::Flowpipe, arr::AbstractVector) = tend(view(array(F), arr))
tspan(F::Flowpipe, arr::UnitRange) = tspan(view(array(F), arr), contiguous=true)
tspan(F::Flowpipe, arr::AbstractVector) = tspan(view(array(F), arr))

# further setops
LazySets.is_intersection_empty(F::Flowpipe{N, <:AbstractLazyReachSet}, Y::LazySet) where {N} = all(X -> _is_intersection_empty(X, Y), array(F))
Base.:⊆(F::Flowpipe, X::LazySet) = all(R ⊆ X for R in F)
Base.:⊆(F::Flowpipe, Y::AbstractLazyReachSet) = all(R ⊆ set(Y) for R in F)

# lazy projection of a flowpipe
function Projection(F::Flowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(F)
    out = map(X -> Projection(X, vars), Xk)
    return Flowpipe(out)
end
Projection(F::Flowpipe; vars) = Projection(F, Tuple(vars))
Projection(F::Flowpipe, vars::AbstractVector{M}) where {M<:Integer} = Projection(F, Tuple(vars))

# membership test
function ∈(x::AbstractVector{N}, fp::Flowpipe{N, <:AbstractLazyReachSet{N}}) where {N}
    return any(R -> x ∈ set(R), array(fp))
end

function ∈(x::AbstractVector{N}, fp::VT) where {N, RT<:AbstractLazyReachSet{N}, VT<:AbstractVector{RT}}
    return any(R -> x ∈ set(R), fp)
end

# =======================================
# Flowpipe composition with a time-shift
# =======================================

"""
    ShiftedFlowpipe{FT<:AbstractFlowpipe, NT<:Number} <: AbstractFlowpipe

Type that lazily represents a flowpipe that has been shifted in time.

### Fields

- `F`  -- original flowpipe
- `t0` -- time shift

### Notes

This type can wrap any concrete subtype of `AbstractFlowpipe`, and the extra
field `t0` is such that the time spans of each reach-set in `F` are shifted
by the amount `t0` (which should be a subtype of `Number`).

A convenience constructor alias `Shift` is given.
"""
struct ShiftedFlowpipe{FT<:AbstractFlowpipe, NT<:Number} <: AbstractFlowpipe
    F::FT
    t0::NT
end

function ShiftedFlowpipe(vec::AbstractVector{<:AbstractLazyReachSet}, t0::Number)
    return ShiftedFlowpipe(Flowpipe(vec), t0)
end

# getter functions
@inline array(fp::ShiftedFlowpipe) = array(fp.F)
@inline flowpipe(fp::ShiftedFlowpipe) = fp.F
@inline time_shift(fp::ShiftedFlowpipe) = fp.t0

# time domain interface
@inline tstart(fp::ShiftedFlowpipe) = tstart(first(fp)) + fp.t0
@inline tend(fp::ShiftedFlowpipe) = tend(last(fp)) + fp.t0
@inline tspan(fp::ShiftedFlowpipe) = TimeInterval(tstart(fp), tend(fp))

@inline tstart(fp::ShiftedFlowpipe, i::Int) = tstart(fp.F[i]) + fp.t0
@inline tend(fp::ShiftedFlowpipe, i::Int) = tend(fp.F[i]) + fp.t0
@inline tspan(fp::ShiftedFlowpipe, i::Int) = TimeInterval(tstart(fp, i), tend(fp, i))

# TODO use interface?
project(fp::ShiftedFlowpipe, vars::AbstractVector) = project(fp, Tuple(vars))
project(fp::ShiftedFlowpipe; vars) = project(fp, Tuple(vars))

function project(fp::ShiftedFlowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(fp)
    # TODO: use projection of the reachsets
    if 0 ∈ vars # projection includes "time"
        # we shift the vars indices by one as we take the Cartesian prod with the time spans
        aux = vars .+ 1
        t0 = time_shift(fp)
        return map(X -> _project(convert(Interval, tspan(X) + t0) × set(X), aux), Xk)
    else
        return map(X -> _project(set(X), vars), Xk)
    end
end

# this method is analogue to project(::AbstractLazyReachSet, vars; check_vars=true)
# TODO add check_vars ?
function project(fp::ShiftedFlowpipe, i::Int, vars::NTuple{D, M}) where {D, M<:Integer}
    t0 = time_shift(fp)
    R = fp[i]
    if 0 ∈ vars
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = vars .+ 1

        Δt = convert(Interval, tspan(R) + t0)
        proj =  _project(Δt × set(R), aux)
    else
        proj = _project(set(R), vars)
    end

    return SparseReachSet(proj, tspan(R) + t0, vars)
end

# TODO: improve using mutable ShiftedFlowpipe struct (so that we can modify t0->t0+t1)
function shift(fp::ShiftedFlowpipe, t1::Number)
    return ShiftedFlowpipe(shift(fp.F, t1), fp.t0)
end

# =====================================
# Flowpipe composition with a lazy map
# =====================================

"""
    MappedFlowpipe{FT<:AbstractFlowpipe, ST} <: AbstractFlowpipe

### Fields

- `F`    -- flowpipe
- `func` -- function representing the map
"""
struct MappedFlowpipe{FT<:AbstractFlowpipe, ST} <: AbstractFlowpipe
    F::FT
    func::ST
end

"""
    Projection(fp::AbstractFlowpipe, vars::NTuple{D, T}) where {D, T<:Integer}

Return the lazy projection of a flowpipe.

### Input

### Output

### Notes

The projection is lazy, and consists of mapping each set
`X` in the flowpipe to `MX`, where `M` is the projection matrix associated with
the given variables `vars`.
"""
function LazySets.Projection(fp::AbstractFlowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    # TODO: assert that vars belong to the variables of the flowpipe
    M = projection_matrix(collect(vars), dim(F), Float64)
    func = @map(x -> M*x)
    return MappedFlowpipe(fp, func)
end

function overapproximate(fp::Flowpipe, args...)
    return Flowpipe(map(R -> overapproximate(R, args...), fp), fp.ext)
end

function overapproximate(fp::VRT, args...) where {RT, VRT<:AbstractVector{RT}}
    return Flowpipe(map(R -> overapproximate(R, args...), fp))
end

# ================================
# Hybrid flowpipe
# ================================

"""
    HybridFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same type, such that they are
contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

The evaluation functions (in time) for this type do not assume that the flowpipes are contiguous in time.
That is, the final time of the `i`-th flowpipe does not match the start time of the `i+1`-th flowpipe.
"""
struct HybridFlowpipe{N, RT<:AbstractReachSet{N}, FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::VectorOfArray{RT, 2, Vector{FT}}
    ext::Dict{Symbol, Any}
end

function HybridFlowpipe(Fk::Vector{FT}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    ext = Dict{Symbol, Any}()
    return HybridFlowpipe{N, RT, FT}(voa, ext)
end

function HybridFlowpipe(Fk::Vector{FT}, ext::Dict{Symbol, Any}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    return HybridFlowpipe{N, RT, FT}(voa, ext)
end

function HybridFlowpipe(Fk::Vector{SFT}) where {N, RT, FT<:Flowpipe{N, RT}, NT<:Number, SFT<:ShiftedFlowpipe{FT, NT}}
    voa = VectorOfArray{RT, 2, Vector{SFT}}(Fk)
    ext = Dict{Symbol, Any}()
    return HybridFlowpipe{N, RT, SFT}(voa, ext)
end

# interface functions
array(fp::HybridFlowpipe) = fp.Fk
flowpipe(fp::HybridFlowpipe) = fp
numtype(::HybridFlowpipe{N}) where {N} = N
setrep(::Type{HybridFlowpipe{N, RT, FT}}) where {N, RT, FT} = setrep(RT)
setrep(::HybridFlowpipe{N, RT, FT}) where {N, RT, FT} = setrep(RT)
rsetrep(::Type{HybridFlowpipe{N, RT, FT}}) where {N, RT, FT} = RT
numrsets(fp::HybridFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::HybridFlowpipe, I::Int...) = getindex(fp.Fk, I...)

function tspan(fp::HybridFlowpipe)
    ti = minimum(tstart, fp)
    tf = maximum(tend, fp)
    return TimeInterval(ti, tf)
end

function Base.similar(fp::HybridFlowpipe{N, RT, FT}) where {N, RT, FT}
    return HybridFlowpipe(Vector{FT}())
end

function overapproximate(fp::HybridFlowpipe, args...)
    return HybridFlowpipe([overapproximate(F, args...) for F in fp], fp.ext)
end

function project(fp::HybridFlowpipe, args...)
    return [project(F, args...) for F in fp]
end

# LazySets interface
function LazySets.ρ(d::AbstractVector, fp::HybridFlowpipe)
    return maximum(ρ(d, F) for F in array(fp))
end

#function LazySets.σ(d::AbstractVector, fp::HybridFlowpipe)
#    error("not implemented")
#end

Base.:⊆(F::HybridFlowpipe, X::LazySet) = all(fp ⊆ X for fp in F)
Base.:⊆(F::HybridFlowpipe, Y::AbstractLazyReachSet) = all(fp ⊆ set(Y) for fp in F)

# evaluation for scalars
function (fp::HybridFlowpipe{N, RT})(t::Number) where {N, RT<:AbstractReachSet{N}}
    if t ∉ tspan(fp)
        throw(ArgumentError("time $t does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        if t ∈ tspan(fk)
            for Ri in fk  # loop over reach-sets for this flowpipe
                if t ∈ tspan(Ri)
                    push!(vec, Ri)
                end
            end
        end
    end
    return vec
end

# evaluation for time intervals
function (fp::HybridFlowpipe{N, RT})(dt::TimeInterval) where {N, RT<:AbstractReachSet{N}}
    if !(dt ⊆ tspan(fp)) # TODO IntervalArithmetic#409
        throw(ArgumentError("time interval $dt does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        if !isdisjoint(dt, tspan(fk))
            for R in fk(dt)  # loop over reach-sets for this flowpipe
                push!(vec, R)
            end
        end
    end
    return vec
end

# ============================================
# Flowpipes not contiguous in time
# ============================================

"""
    MixedFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same time, such that they are
not necessarily contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

This type does not assume that the flowpipes are contiguous in time.
"""
struct MixedFlowpipe{N, RT<:AbstractReachSet{N}, FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::VectorOfArray{RT, 2, Vector{FT}}
    ext::Dict{Symbol, Any}
end

function MixedFlowpipe(Fk::Vector{FT}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    ext = Dict{Symbol, Any}()
    return MixedFlowpipe{N, RT, FT}(voa, ext)
end

function MixedFlowpipe(Fk::Vector{FT}, ext::Dict{Symbol, Any}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    return MixedFlowpipe{N, RT, FT}(voa, ext)
end

# interface functions
array(fp::MixedFlowpipe) = fp.Fk
flowpipe(fp::MixedFlowpipe) = fp
setrep(::Type{MixedFlowpipe{N, RT, FT}}) where {N, RT, FT} = RT
numrsets(fp::MixedFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::MixedFlowpipe, I::Int...) = getindex(fp.Fk, I...)

@inline tspan(fp::MixedFlowpipe) = error("for mixed flowpipes you should specify the index, as in `tspan(fp, i)`")
@inline tspan(fp::MixedFlowpipe, i::Int) = tspan(fp[i])

function Base.similar(fp::MixedFlowpipe{N, RT, FT}) where {N, RT, FT}
    return MixedFlowpipe(Vector{FT}())
end

function (fp::MixedFlowpipe)(t::Number)
    error("not implemented yet")
end

function (fp::MixedFlowpipe)(dt::TimeInterval)
    error("not implemented yet")
end

function overapproximate(fp::MixedFlowpipe, args...)
    return MixedFlowpipe([overapproximate(Fi, args...) for Fi in fp], fp.ext)
end

function project(fp::MixedFlowpipe, args...)
    return [project(F, args...) for F in fp]
end

# LazySets interface
function LazySets.ρ(d::AbstractVector, fp::MixedFlowpipe)
    return maximum(ρ(d, F) for F in array(fp))
end

function LazySets.σ(d::AbstractVector, fp::MixedFlowpipe)
    error("not implemented")
end

# ============================================
# Hybrid flowpipe of possibly different types
# ============================================

"""
    PartitionedFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of possibly different types.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes
"""
struct PartitionedFlowpipe{T, S<:Tuple} <: AbstractFlowpipe # TODO: ask <:AbstractFlowpipe for each element in the tuple..?
    Fk::ArrayPartition{T, S}
    ext::Dict{Symbol, Any}
end
