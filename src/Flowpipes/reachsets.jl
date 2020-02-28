export ReachSet,
       SparseReachSet

export project,
       set,
       setrep,
       tstart,
       tend,
       tspan

# ================================
# Reach set interfaces
# ================================

"""
    AbstractReachSet{N}

Abstract type for all reach-sets types.

### Notes

A reach-set is a set representation `X` associated to a given time span `Δt`.

In its simplest form, we represent reach-sets with a struct that wraps the tuple
`(X, Δt)`, where `X` corresponds to a geometric set, eg. a polytope, and `Δt` is
the interval with the time span associated to this reach-set. Concrete subtypes
of `AbstractReachSet` may represent the set `X` in different ways, or add extra
information associated to variables, etc.

This type is parametric in `N`. This parameter is used in some subtypes for the
type of numerical coefficient used in the representation (typically, `N = Float64`).

Each concrete subtype should implement the following methods:

- `set`    -- return the geometric set
- `setrep` -- return the type of the set representation
- `tspan`  -- return the time interval span
- `tstart` -- return the initial time
- `tend`   -- return the final time
- `dim`    -- ambient dimension of the reach-set
"""
abstract type AbstractReachSet{N} end

"""
    set(R::AbstractReachSet)

Return the geometric set represented by this reach-set.

## Input

- `R` -- reach-set

## Output

The set wrapped by the given reach-set.
"""
function set(::AbstractReachSet) end

"""
    setrep(R::AbstractReachSet)

Return the type of the set representation of this reach-set.

## Input

- `R` -- reach-set

## Output

Type of the set representation of the given reach-set.
"""
setrep(::AbstractReachSet)

"""
    tspan(R::AbstractReachSet)

Return time span of this reach-set.

## Input

- `R` -- reach-set

## Output

The interval representing the time span of the given reach-set.
"""
function tspan(::AbstractReachSet) end

"""
    tstart(R::AbstractReachSet)

Return the initial time of this reach-set.

## Input

- `R` -- reach-set

## Output

A float representing the initial time of the given reach-set.
"""
function tstart(::AbstractReachSet) end

"""
    tend(R::AbstractReachSet)

Return the final time of this reach-set.

## Input

- `R` -- reach-set

## Output

A float representing the final time of the given reach-set.
"""
function tend(::AbstractReachSet) end

"""
    dim(R::AbstractReachSet)

Return the ambient dimension of the reach-set.

## Input

- `R` -- reach-set

## Output

An integer corresponding to the ambient dimension of the given reach-set.
"""
function LazySets.dim(::AbstractReachSet) end

"""
    AbstractLazyReachSet{N} <: AbstractReachSet{N}

Abstract type for all reach-set types that use a `LazySet` for the underlying
set representation.

### Notes

An `AbstractLazyReachSet` is the interface for reach-sets such that the geometric
set is represented by any subtype of `LazySet`.

This types implements the `LazySets` interface, namely support function (`ρ`),
support vector (`σ`) and ambient dimension (`dim`) functions. Hence, these
functions directly apply to concrete subtypes of an `AbstractLazyReachSet`. The
set wrapped by this type is obtained through `set(R)`.

In addition to the functions inherited from `AbstractReachSet`, the following
are available:

- `project` -- projection of a reach-set along the given variables in `vars`
- `vars`    -- tuple of integers associated to the variables of the given reach-set
"""
abstract type AbstractLazyReachSet{N} <: AbstractReachSet{N} end

# Implement LazySets interface
LazySets.ρ(d::AbstractVector, R::AbstractLazyReachSet) = ρ(d, set(R))
LazySets.σ(d::AbstractVector, R::AbstractLazyReachSet) = σ(d, set(R))
LazySets.dim(R::AbstractLazyReachSet) = dim(set(R))

# Expose common LazySets operations
# radius, etc TODO: #24


"""
    project(R::AbstractLazyReachSet, vars::NTuple{D, Int}) where {D}

Projects a reach-set onto the subspace spanned by the given variables.

### Input

- `R`    -- reach-set
- `vars` -- tuple of variables for the projection

### Notes

This function can be used to project a reach-set onto a lower-dimensional
sub-space. The projection is lazy, and consists of mapping `X = set(R)` to `MX`,
where `M` is the projection matrix associated with the given variables `vars`.
"""
function project(R::AbstractLazyReachSet, vars::NTuple{D, Int}; lazy::Bool=true) where {D}
    if !(vars ⊆ vars(R))
        throw(ArgumentError("the variables `vars` do not belong to the variables " *
                " of this reach-set, $(vars(R))"))
    end
    if lazy
        X = set(R)
        return project(X, vars, LinearMap)
    else
        error("the concrete projection is not implemented yet")
    end
end

# handle generic vars vector
function project(R::AbstractLazyReachSet, vars::AbstractVector; lazy::Bool=true)
    vars = Tuple(vars) # Tuple(vi for vi in vcat(vars...))
    return project(R, vars; lazy=lazy)
end

"""
    AbstractTaylorModelReachSet{N}

Abstract type for all reach sets types that represent a Taylor model.

### Notes

The parameter `N` refers to the numerical type of the representation.
"""
abstract type AbstractTaylorModelReachSet{N} <: AbstractReachSet{N} end

# ================================
# Reach set
# ================================

"""
    ReachSet{N, ST<:LazySet{N}} <: AbstractLazyReachSet{N}

Type that wraps a reach-set using a `LazySet` as underlying representation.

### Fields

- `X`  -- set
- `Δt` -- time interval

### Notes

A `ReachSet` is a struct representing (an approximation of) the reachable states
for a given time interval. The type of the representation is `ST`, which may be
any subtype LazySet. For efficiency reasons, `ST` should be concretely typed.

By assumption the coordinates in this reach-set are associated to the integers
`1, …, n`.
"""
struct ReachSet{N, ST<:LazySet{N}} <: AbstractLazyReachSet{N}
    X::ST
    Δt::IA.Interval{Float64}
end

# abstract reach set interface functions
set(R::ReachSet) = R.X
setrep(R::ReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
tstart(R::ReachSet) = inf(R.Δt)
tend(R::ReachSet) = sup(R.Δt)
tspan(R::ReachSet) = R.Δt
LazySets.dim(R::ReachSet) = dim(R.X)
vars(R::ReachSet) = Tuple(Base.OneTo(dim(R.X)),)

# ================================
# Sparse reach set
# ================================

"""
    SparseReachSet{N, ST<:LazySet{N}, D} <: AbstractReachSet{N}

Type that wraps a reach-set using a `LazySet` as underlying representation,
together with a tuple of variables associated to this reach-set.

### Fields

- `X`    -- set
- `Δt`   -- time interval
- `vars` -- variables represented by the set `X`

### Notes

A `SparseReachSet` is a struct representing (an approximation of) the reachable
states for a given time interval. The type of the representation is `ST`, which
may be any subtype of `LazySet` (ideally, concrete). Moreover, this type also
stores information about the variables (also named coordinates, or by abuse of
notation, *dimensions*) corresponding to the set `X`.

For instance in the ambient space `n=5`, one may have a `SparseReachSet` whose
variables tuple is `vars = (4, 5, 6)`, i.e. representing a three-dimensional
projection of the full-dimensional reach-set. In consequence, the dimension of
`X` doesn't match the length of `vars`, in general

In this type, the parameter `N` represents the numerical type of the `LazySet`
(typically, `Float64`), the type `ST` represents the set representation used,
and `D` denotes the dimension of this sparse reach set. Note that, in contrast
to `ReachSet`, for `SparseReachSet` the number of dimensions is part of the type
information.
"""
struct SparseReachSet{N, ST<:LazySet{N}, D} <: AbstractLazyReachSet{N}
    X::ST
    Δt::IA.Interval{Float64}
    vars::NTuple{D, Int}
end

# interface functions
set(R::SparseReachSet) = R.X
setrep(R::SparseReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
tstart(R::SparseReachSet) = inf(R.Δt)
tend(R::SparseReachSet) = sup(R.Δt)
tspan(R::SparseReachSet) = R.Δt
LazySets.dim(R::SparseReachSet{N, ST, D}) where {N, ST<:LazySet{N}, D} = D
vars(R::SparseReachSet) = R.vars

# constructor from vector of dimensions
function SparseReachSet(X::ST, Δt::IA.Interval{Float64},
                        vars::AbstractVector) where {N, ST<:LazySet{N}}
    SparseReachSet(X, Δt, Tuple(vars))
end

# ================================
# Taylor model reach set
# ================================

#=
struct TaylorModelReachSet{N, ST, D} <: AbstractTaylorModelReachSet{N}
    X::ST
    Δt::IA.Interval{Float64}
end

...
=#
