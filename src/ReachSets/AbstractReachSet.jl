# ================================================================
# Reach set interfaces
# ================================================================

"""
    AbstractReachSet{N}

Abstract type for all reach-sets types.

### Notes

A reach-set is a set representation `X` associated to a given time span `Δt`.

In its simplest form, we represent reach-sets with a struct that wraps the tuple
`(X, Δt)`, where `X` corresponds to a geometric set, eg. a polytope, and `Δt` is
the interval with the time span associated to this reach-set.

This type is parametric in `N`. The parameter `N` represents for the type of numerical
coefficient used in the representation (typically, `N = Float64`).

Although concrete subtypes of `AbstractReachSet` may represent the set `X` in
different ways, or carry additional information as an extra type field, they
should all implement the `AbstractReachSet` interface to enable shared functionality
for reach-set types. In particular, each concrete subtype should implement the
following methods:

- `set`    -- return the geometric set
- `setrep` -- return the type of the set representation
- `tspan`  -- return the time interval span
- `tstart` -- return the initial time
- `tend`   -- return the final time
- `dim`    -- return ambient dimension of the reach-set
"""
abstract type AbstractReachSet{N} end

# convenience union for dispatch on structs that behave like a set
const SetOrReachSet = Union{LazySet,UnionSet,UnionSetArray,IA.Interval,IA.IntervalBox,
                            AbstractReachSet}

"""
    basetype(T::Type{<:AbstractReachSet})

Return the base type of the given reach-set type (i.e., without type parameters).

### Input

- `T` -- reach-set type, used for dispatch

### Output

The base type of `T`.
"""
basetype(T::Type{<:AbstractReachSet}) = Base.typename(T).wrapper

"""
    set(R::AbstractReachSet)

Return the geometric set represented by this reach-set.

### Input

- `R` -- reach-set

### Output

The set wrapped by the given reach-set.
"""
function set(::AbstractReachSet) end

# no-op
set(X::LazySet) = X

"""
    setrep(R::AbstractReachSet)

Return the type of the set representation of this reach-set.

### Input

- `R` -- reach-set

### Output

Type of the set representation of the given reach-set.
"""
setrep(::AbstractReachSet)

"""
    tstart(R::AbstractReachSet)

Return the initial time of this reach-set.

### Input

- `R` -- reach-set

### Output

A float representing the initial time of the given reach-set.
"""
function tstart(::AbstractReachSet) end

"""
    tend(R::AbstractReachSet)

Return the final time of this reach-set.

### Input

- `R` -- reach-set

### Output

A float representing the final time of the given reach-set.
"""
function tend(::AbstractReachSet) end

"""
    tspan(R::AbstractReachSet)

Return time span of this reach-set.

### Input

- `R` -- reach-set

### Output

The interval representing the time span of the given reach-set.
"""
function tspan(::AbstractReachSet) end

# if `contiguous` is true, we assume that the time span of each reach-set in this slice
# is sorted, e.g. F[1:4] so it suffices to look at the first and last reach-sets
function tstart(F::VRT; contiguous=false) where {RT<:AbstractReachSet,VRT<:AbstractVector{RT}}
    if contiguous
        t0 = tstart(first(F))
    else
        t0 = minimum(tstart, F)
    end
    return t0
end

# if `contiguous` is true, we assume that the time span of each reach-set in this slice
# is sorted, e.g. F[1:4] so it suffices to look at the first and last reach-sets
function tend(F::VRT; contiguous=false) where {RT<:AbstractReachSet,VRT<:AbstractVector{RT}}
    if contiguous
        tf = tend(last(F))
    else
        tf = maximum(tend, F)
    end
    return tf
end

# if `contiguous` is true, we assume that the time span of each reach-set in this slice
# is sorted, e.g. F[1:4] so it suffices to look at the first and last reach-sets
function tspan(F::VRT; contiguous=false) where {RT<:AbstractReachSet,VRT<:AbstractVector{RT}}
    if contiguous
        t0 = tstart(first(F))
        tf = tend(last(F))
    else
        t0 = minimum(tstart, F)
        tf = maximum(tend, F)
    end
    return TimeInterval(t0, tf)
end

"""
    dim(R::AbstractReachSet)

Return the ambient dimension of the reach-set.

### Input

- `R` -- reach-set

### Output

An integer corresponding to the ambient dimension of the given reach-set.
"""
function LazySets.dim(::AbstractReachSet) end

# for internal use (dimensional checks in solve)
_dim(X::AbstractReachSet) = dim(X)

"""
    copy(R::AbstractReachSet)

Return a copy of the given reach-set.

### Input

- `R`  -- reach-set

### Output

A new reach-set of the sam type and the same field values as `R`.
"""
Base.copy(R::AbstractReachSet) = deepcopy(R)

"""
    shift(R::AbstractReachSet, t0::Number)

Perform a time-shift of the given reach-set.

### Input

- `R`  -- reach-set
- `t0` -- number that corresponds to the time-shift

### Output

A new reach-set of the same type of `R` such that its time-span has been shifted
by `t0`.
"""
function shift(R::AbstractReachSet, t0::Number) end

# ------------------------------
# Methods to check disjointness
# ------------------------------

# fallback uses  internal function _is_intersection_empty, which admit a pre-processing
# step for the reach-set / algorithm choice
@commutative function isdisjoint(R::AbstractReachSet, Y::SetOrReachSet,
                                 method::AbstractDisjointnessMethod=FallbackDisjointness())
    return _is_intersection_empty(R, Y, method)
end

# disambiguations
function isdisjoint(R1::AbstractReachSet, R2::AbstractReachSet,
                    method::AbstractDisjointnessMethod=FallbackDisjointness())
    return _is_intersection_empty(R1, R2, method)
end

# vector of reach-sets
@commutative function isdisjoint(R::AbstractVector{<:AbstractReachSet}, Y::SetOrReachSet,
                                 method::AbstractDisjointnessMethod=FallbackDisjointness())
    return all(Ri -> _is_intersection_empty(Ri, Y, method), R)
end

# internal implementation
function _is_intersection_empty(R::AbstractReachSet, Y::AbstractReachSet,
                                method::FallbackDisjointness)
    return isdisjoint(set(R), set(Y))
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet,
                                             method::FallbackDisjointness)
    return isdisjoint(set(R), Y)
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet,
                                             method::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return isdisjoint(set(Z), Y)
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, method::BoxEnclosure)
    H = overapproximate(R, Hyperrectangle)
    return isdisjoint(set(H), Y)
end

# used for disjointness check in continuous post-operators
function _is_intersection_empty(P::HPolyhedron, Z::Zonotope, ::BoxEnclosure)
    return isdisjoint(P, box_approximation(Z))
end

# in this method we assume that the intersection is non-empty
@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, method::Dummy)
    return false
end

@commutative function _is_intersection_empty(R::AbstractReachSet,
                                             Y::UnionSet{N,<:Interval{N},<:Interval{N}},
                                             method::FallbackDisjointness) where {N}
    return isdisjoint(set(R), Y)
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::UnionSetArray{N,<:Interval{N}},
                                             method::FallbackDisjointness) where {N}
    return isdisjoint(set(R), Y)
end

# ------------------------------
# Concrete intersection methods
# ------------------------------

function intersection(R::AbstractReachSet, S::AbstractReachSet)
    return _intersection(set(R), set(S), FallbackIntersection())
end
function _intersection(R::AbstractReachSet, S::AbstractReachSet, method::AbstractIntersectionMethod)
    return _intersection(set(R), set(S), method)
end

# fallback methods for reach-sets
@commutative function _intersection(R::AbstractReachSet, X::LazySet,
                                    method::AbstractIntersectionMethod)
    return _intersection(set(R), X, method)
end

# ------------------------------
# Methods to overapproximate
# ------------------------------

box_approximation(R::AbstractReachSet; kwargs...) = overapproximate(R, Hyperrectangle; kwargs...)
