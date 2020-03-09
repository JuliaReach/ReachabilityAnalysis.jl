# types
export ReachSet,
       SparseReachSet,
       Projection

# methods
export set,
       setrep,
       tstart,
       tend,
       tspan,
       project,
       shift

# method extensions
import LazySets: dim

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

# ============================================================
# AbstractLazyReachSet: reach-set that behaves like a LazySet
# ============================================================

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

The following functions should be implemented by any concrete subtype:

- `reconstruct` -- create a new instance of the given reach-set with a different
                   set representation but sharing the other fields, i.e. the same
                   time span (and the same for other fields, if applicable)

In addition to the functions inherited from `AbstractReachSet`, the following
are available:

- `project`  -- projection of a reach-set
- `shift`    -- time-shift of a reach-set
- `vars_idx` -- tuple of integers associated to the variables of the given reach-set
"""
abstract type AbstractLazyReachSet{N} <: AbstractReachSet{N} end

# Implement LazySets interface
LazySets.ρ(d::AbstractVector, R::AbstractLazyReachSet) = ρ(d, set(R))
LazySets.σ(d::AbstractVector, R::AbstractLazyReachSet) = σ(d, set(R))
LazySets.dim(R::AbstractLazyReachSet) = dim(set(R))

# Expose common LazySets operations
LazySets.constraints_list(R::AbstractLazyReachSet) = constraints_list(set(R))
LazySets.vertices_list(R::AbstractLazyReachSet) = vertices_list(set(R))
function LazySets.LinearMap(M::Union{AbstractMatrix, Number}, R::AbstractLazyReachSet)
    return reconstruct(R, LinearMap(M, set(R)))
end
function LazySets.linear_map(M::AbstractMatrix, R::AbstractLazyReachSet)
    return reconstruct(R, linear_map(M, set(R)))
end

# handle generic vars vector
function project(R::AbstractLazyReachSet, vars::AbstractVector{M}) where {M<:Integer}
    return project(R, Tuple(vars))
end

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
`1, …, n`. The function `vars_idx` returns such tuple.
"""
struct ReachSet{N, ST<:LazySet{N}} <: AbstractLazyReachSet{N}
    X::ST
    Δt::IA.Interval{Float64}
end

# abstract reach set interface functions
set(R::ReachSet) = R.X
setrep(R::ReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
setrep(::Type{ReachSet{N, ST}}) where {N, ST<:LazySet{N}} = ST
tstart(R::ReachSet) = inf(R.Δt)
tend(R::ReachSet) = sup(R.Δt)
tspan(R::ReachSet) = R.Δt
LazySets.dim(R::ReachSet) = dim(R.X)
vars_idx(R::ReachSet) = Tuple(Base.OneTo(dim(R.X)),)

# no-op
function Base.convert(T::Type{ReachSet{N, ST}}, R::ReachSet{N, ST}) where {N, ST<:LazySet{N}}
    return R
end

# convert to a reach-set with another underlying set representation
function Base.convert(T::Type{ReachSet{N, LT}}, R::ReachSet{N, ST}) where {N, ST<:LazySet{N}, LT<:LazySet{N}}
    X = convert(LT, set(R))
    return ReachSet(X, tspan(R))
end

function shift(R::ReachSet, t0::Number)
    return ReachSet(set(R), tspan(R) + t0)
end

function reconstruct(R::ReachSet, Y::LazySet)
    return ReachSet(Y, tspan(R))
end

# ================================
# Sparse reach set
# ================================

"""
    SparseReachSet{N, ST<:LazySet{N}, D} <: AbstractReachSet{N}

Type that wraps a reach-set using a `LazySet` as underlying representation,
together with a tuple of variables associated to this reach-set.

### Fields

- `X`        -- set
- `Δt`       -- time interval
- `vars_idx` -- tuple of variable indices represented by the set `X`

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
    # TODO: inner constructor that the dimension of vars matches that of X ?
end

# interface functions
set(R::SparseReachSet) = R.X
setrep(R::SparseReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
tstart(R::SparseReachSet) = inf(R.Δt)
tend(R::SparseReachSet) = sup(R.Δt)
tspan(R::SparseReachSet) = R.Δt
LazySets.dim(R::SparseReachSet{N, ST, D}) where {N, ST<:LazySet{N}, D} = D
vars_idx(R::SparseReachSet) = R.vars

# constructor from vector of dimensions
function SparseReachSet(X::ST, Δt::IA.Interval{Float64},
                        vars::AbstractVector) where {N, ST<:LazySet{N}}
    SparseReachSet(X, Δt, Tuple(vars))
end

function shift(R::SparseReachSet, t0::Number)
    return SparseReachSet(set(R), tspan(R) + t0, vars_idx(R))
end

function reconstruct(R::SparseReachSet, Y::LazySet)
    return SparseReachSet(Y, tspan(R), vars_idx(R))
end

"""
    project(R::Union{ReachSet, SparseReachSet}, vars::NTuple{D, Int};
            [check_vars]::Bool=true) where {D}

Projects a reach-set onto the subspace spanned by the given variables.

### Input

- `R`          -- reach-set
- `vars`       -- tuple of variables for the projection
- `check_vars` -- (optional, default: `true`) if `true`, check that the given variable
                  indices `vars` are a subset of the variables of `R`

### Output

A `SparseReachSet` whose variable indices are given by `vars`.

The type of the new reach-set depends on the type of the reach-set `R`:

- If `R` contains a hyperrectangular set, the output is a hyperrectangle.
- If `R` contains a zonotopic set, the output is a zonotope.
- Otherwise, the return type is a polytope either in constraint representation
  or in vertex representation, depending on the dimension and the properties of `M`.
  For details, see `LazySets.project`.

### Notes

This function can be used to project a reach-set onto a lower-dimensional
sub-space. The projection is concrete, and it consists of mapping the reach-set
`X = set(R)` to a new reach-set through to `MX`, where `M` is the projection matrix
associated with the given variables `vars`.

To project onto the time variable, use the index `0`. For instance, `(0, 1)` projects
onto the time variable and the first variable in `R`.
"""
function project(R::AbstractLazyReachSet, vars::NTuple{D, M};
                 check_vars::Bool=true) where {D, M<:Integer}

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(vars, 0) ⊆ vars_idx(R))
        throw(ArgumentError("the variables $vars do not belong to the variables " *
                            " of this reach-set, $(vars_idx(R))"))
    end

    if 0 ∈ vars
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = vars .+ 1
        Δt = convert(Interval, tspan(R))
        proj =  _project(Δt × set(R), aux)
    else
        proj = _project(set(R), vars)
    end

    return SparseReachSet(proj, tspan(R), vars)
end

# lazy projection of a reach-set
function Projection(R::AbstractLazyReachSet, vars::NTuple{D, M},
                    check_vars::Bool=true) where {D, M<:Integer}

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(vars, 0) ⊆ vars_idx(R))
        throw(ArgumentError("the variables $vars do not belong to the variables " *
                            " of this reach-set, $(vars_idx(R))"))
    end

    if 0 ∈ vars
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = vars .+ 1
        Δt = convert(Interval, tspan(R))
        proj =  _Projection(Δt × set(R), aux)
    else
        proj = _Projection(set(R), vars)
    end

    return SparseReachSet(proj, tspan(R), vars)
end

# ================================
# Taylor model reach set
# ================================

using TaylorModels: TaylorModel1, TaylorN

"""
    AbstractTaylorModelReachSet{N}

Abstract type for all reach sets types that represent a Taylor model.

### Notes

The parameter `N` refers to the numerical type of the representation.
"""
abstract type AbstractTaylorModelReachSet{N} <: AbstractReachSet{N} end

"""
    TaylorModelReachSet{N} <: AbstractTaylorModelReachSet{N}

Taylor model reach-set represented as a vector taylor models in one variable
(namely, the "time" variable) whose coefficients are multivariate polynomials
(namely in the "space" variables).

### Notes

The parameter `N` refers to the numerical type of the representation.
"""
struct TaylorModelReachSet{N, D<:Integer} <: AbstractTaylorModelReachSet{N}
    X::Vector{TaylorModel1{TaylorN{N}, N}}
    Δt::TimeInterval
    #dim::D
end

# interface functions
set(R::TaylorModelReachSet) = R.X
setrep(R::TaylorModelReachSet{N}) where {N} = Vector{TaylorModel1{TaylorN{N}, N}}
tstart(R::TaylorModelReachSet) = inf(R.Δt)
tend(R::TaylorModelReachSet) = sup(R.Δt)
tspan(R::TaylorModelReachSet) = R.Δt
LazySets.dim(R::TaylorModelReachSet{N, D}) where {N, D} = R.dim # TODO (!)
vars_idx(R::TaylorModelReachSet) = Tuple(Base.OneTo(dim(R.X)),)

function shift(R::TaylorModelReachSet, t0::Number)
    return TaylorModelReachSet(set(R), tspan(R) + t0)
end
