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
basetype(fp::AbstractFlowpipe) = basetype(typeof(fp))

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

# fallback implementations of set getter functions

"""
    set(fp::AbstractFlowpipe)

Return the geometric set represented by this flowpipe as the union of reach-sets.

### Input

- `fp`  -- flowpipe

### Output

The set union of the array of reach-sets of the flowpipe.

## Notes

To retrieve the array of sets stored in the flowpipe use `array(fp)`. To get
a set at a particular index, use `set(F[ind])` or `set(F, ind)`.
"""
set(fp::AbstractFlowpipe) = UnionSetArray([set(R) for R in array(fp)])

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

"""
    set(fp::AbstractFlowpipe, ind::AbstractVector)

Return the union of set represented by this flowpipe at the given indices.

## Input

- `fp`  -- flowpipe
- `ind` -- vector of indices

## Output

The set union stored in the flowpipe at the given indices.
"""
set(fp::AbstractFlowpipe, ind::AbstractVector) = UnionSetArray([set(fp, i) for i in ind])

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

- `fp` -- flowpipe

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

"""
    vars(fp::AbstractFlowpipe)

Return the tuple of variable indices of the flowpipe.

### Input

- `fp` -- flowpipe

### Output

Tuple of integers with the variable indices of the flowpipe, typically ``1, 2, …, n``
where ``n`` is the dimension of the flowpipe.

### Notes

The fallback implementation assumes first reach-set is representative.
"""
vars(fp::AbstractFlowpipe) = vars(first(fp))

# indexing with ranges or with vectors of integers
Base.getindex(fp::AbstractFlowpipe, i::Int) = getindex(array(fp), i)
Base.getindex(fp::AbstractFlowpipe, i::Number) = getindex(array(fp), convert(Int, i))
Base.getindex(fp::AbstractFlowpipe, I::AbstractVector) = getindex(array(fp), I)

# ------------------------------
# Methods to check disjointness
# ------------------------------

# interface
is_intersection_empty(F::AbstractFlowpipe, Y::SetOrReachSet, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> _is_intersection_empty(X, Y, method), array(F))
is_intersection_empty(Y::SetOrReachSet, F::AbstractFlowpipe, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> _is_intersection_empty(X, Y, method), array(F))
Base.isdisjoint(F::AbstractFlowpipe, Y::SetOrReachSet, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> _is_intersection_empty(X, Y, method), array(F))
Base.isdisjoint(Y::SetOrReachSet, F::AbstractFlowpipe, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> _is_intersection_empty(X, Y, method), array(F))

is_intersection_empty(F1::AbstractFlowpipe, F2::AbstractFlowpipe, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> isdisjoint(F1, X, method), array(F2))
Base.isdisjoint(F1::AbstractFlowpipe, F2::AbstractFlowpipe, method::AbstractDisjointnessMethod=FallbackDisjointness()) = all(X -> isdisjoint(F1, X, method), array(F2))

# flowpipe and vector of reach-sets TODO

# internal
_is_intersection_empty(F::AbstractFlowpipe, Y::SetOrReachSet, method::AbstractDisjointnessMethod) = all(X -> _is_intersection_empty(X, Y, method), array(F))

# ------------------------------
# Methods to check inclusion
# ------------------------------

Base.:⊆(F::AbstractFlowpipe, X::LazySet) = all(R ⊆ X for R in F)
Base.:⊆(F::AbstractFlowpipe, Y::AbstractLazyReachSet) = all(R ⊆ set(Y) for R in F)
