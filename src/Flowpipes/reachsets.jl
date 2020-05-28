# method extensions
import LazySets: dim, overapproximate
import Base: ∈

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
function tstart(F::VRT; contiguous=false) where {RT<:AbstractReachSet, VRT<:AbstractVector{RT}}
    if contiguous
        t0 = tstart(first(F))
    else
        t0 = minimum(tstart, F)
    end
    return t0
end

# if `contiguous` is true, we assume that the time span of each reach-set in this slice
# is sorted, e.g. F[1:4] so it suffices to look at the first and last reach-sets
function tend(F::VRT; contiguous=false) where {RT<:AbstractReachSet, VRT<:AbstractVector{RT}}
    if contiguous
        tf = tend(last(F))
    else
        tf = maximum(tend, F)
    end
    return tf
end

# if `contiguous` is true, we assume that the time span of each reach-set in this slice
# is sorted, e.g. F[1:4] so it suffices to look at the first and last reach-sets
function tspan(F::VRT; contiguous=false) where {RT<:AbstractReachSet, VRT<:AbstractVector{RT}}
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
- `vars`     -- tuple of integers associated to the variables of the given reach-set
"""
abstract type AbstractLazyReachSet{N} <: AbstractReachSet{N} end

# Implement LazySets interface
LazySets.ρ(d::AbstractVector, R::AbstractLazyReachSet) = ρ(d, set(R))
LazySets.σ(d::AbstractVector, R::AbstractLazyReachSet) = σ(d, set(R))
LazySets.dim(R::AbstractLazyReachSet) = dim(set(R))

# Expose common LazySets operations
LazySets.constraints_list(R::AbstractLazyReachSet) = constraints_list(set(R))
LazySets.vertices_list(R::AbstractLazyReachSet) = vertices_list(set(R))
Base.:⊆(R::AbstractLazyReachSet, X::LazySet) = ⊆(set(R), X)
Base.:⊆(X::LazySet, R::AbstractLazyReachSet) = ⊆(X, set(R))
Base.:⊆(R::AbstractLazyReachSet, S::AbstractLazyReachSet) = ⊆(set(R), set(S))
LazySets.area(R::AbstractLazyReachSet) = area(set(R))
LazySets.volume(R::AbstractLazyReachSet) = volume(set(R))
Base.convert(::Type{ST}, R::AbstractLazyReachSet) where {ST<:LazySet} = convert(ST, set(R))
Base.convert(::Type{<:IntervalBox}, R::AbstractLazyReachSet) = convert(IntervalBox, set(R))
complement(R::AbstractLazyReachSet) = reconstruct(R, complement(set(R)))

# forward to internal function _is_intersection_empty, which admit a pre-processing
# step for the reach-set / algorithm choice
LazySets.is_intersection_empty(R::AbstractReachSet, Y::LazySet) = _is_intersection_empty(R, Y)

function LazySets.LinearMap(M::Union{AbstractMatrix, Number}, R::AbstractLazyReachSet)
    return reconstruct(R, LinearMap(M, set(R)))
end

function LazySets.linear_map(M::AbstractMatrix, R::AbstractLazyReachSet)
    return reconstruct(R, linear_map(M, set(R)))
end

function LazySets.overapproximate(R::AbstractLazyReachSet, func)
    return reconstruct(R, overapproximate(set(R), func))
end

# handle generic vars vector
function project(R::AbstractLazyReachSet, vars::AbstractVector{M}) where {M<:Integer}
    return project(R, Tuple(vars))
end

# handle generic kwargs vars
function project(R::AbstractLazyReachSet; vars)
    return project(R, Tuple(vars))
end

# membership test
function ∈(x::AbstractVector{N}, R::AbstractLazyReachSet{N}) where {N}
    return ∈(x, set(R))
end

# ================================================================
# Reach set
# ================================================================

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
`1, …, n`. The function `vars` returns such tuple.
"""
struct ReachSet{N, ST<:LazySet{N}} <: AbstractLazyReachSet{N}
    X::ST
    Δt::TimeInterval
end

# abstract reach set interface functions
set(R::ReachSet) = R.X
setrep(R::ReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
setrep(::Type{ReachSet{N, ST}}) where {N, ST<:LazySet{N}} = ST
tstart(R::ReachSet) = inf(R.Δt)
tend(R::ReachSet) = sup(R.Δt)
tspan(R::ReachSet) = R.Δt
dim(R::ReachSet) = dim(R.X)
vars(R::ReachSet) = Tuple(Base.OneTo(dim(R.X)),)

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

# ================================================================
# Sparse reach set
# ================================================================

"""
    SparseReachSet{N, ST<:LazySet{N}, D} <: AbstractReachSet{N}

Type that wraps a reach-set using a `LazySet` as underlying representation,
together with a tuple of variables associated to this reach-set.

### Fields

- `X`        -- set
- `Δt`       -- time interval
- `vars` -- tuple of variable indices represented by the set `X`

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
setrep(::SparseReachSet{N, ST}) where {N, ST<:LazySet{N}} = ST
setrep(::Type{<:SparseReachSet{N, ST}}) where {N, ST<:LazySet{N}} = ST
tstart(R::SparseReachSet) = inf(R.Δt)
tend(R::SparseReachSet) = sup(R.Δt)
tspan(R::SparseReachSet) = R.Δt
dim(R::SparseReachSet{N, ST, D}) where {N, ST<:LazySet{N}, D} = D
vars(R::SparseReachSet) = R.vars

# constructor from vector of dimensions
function SparseReachSet(X::ST, Δt::IA.Interval{Float64},
                        vars::AbstractVector) where {N, ST<:LazySet{N}}
    SparseReachSet(X, Δt, Tuple(vars))
end

function shift(R::SparseReachSet, t0::Number)
    return SparseReachSet(set(R), tspan(R) + t0, vars(R))
end

function reconstruct(R::SparseReachSet, Y::LazySet)
    return SparseReachSet(Y, tspan(R), vars(R))
end

"""
    project(R::Union{ReachSet, SparseReachSet}, variables::NTuple{D, Int};
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
function project(R::AbstractLazyReachSet, variables::NTuple{D, M};
                 check_vars::Bool=true) where {D, M<:Integer}

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(variables, 0) ⊆ vars(R))
        throw(ArgumentError("the variables $vars do not belong to the variables " *
                            " of this reach-set, $(vars(R))"))
    end

    if 0 ∈ variables
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = variables .+ 1
        Δt = convert(Interval, tspan(R))
        proj =  _project(Δt × set(R), aux)
    else
        proj = _project(set(R), variables)
    end

    return SparseReachSet(proj, tspan(R), variables)
end

#=
function project(R::SparseReachSet, variables::NTuple{D, M};
                 check_vars::Bool=true) where {D, M<:Integer}

# TODO
end
=#

# lazy projection of a reach-set
function Projection(R::AbstractLazyReachSet, variables::NTuple{D, M},
                    check_vars::Bool=true) where {D, M<:Integer}

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(variables, 0) ⊆ vars(R))
        throw(ArgumentError("the variables $variables do not belong to the variables " *
                            " of this reach-set, $(vars(R))"))
    end

    if 0 ∈ variables
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = variables .+ 1
        Δt = convert(Interval, tspan(R))
        proj =  _Projection(Δt × set(R), aux)
    else
        proj = _Projection(set(R), variables)
    end

    return SparseReachSet(proj, tspan(R), variables)
end

# handle generic vars vector
function Projection(R::AbstractLazyReachSet, vars::AbstractVector{M}) where {M<:Integer}
    return Projection(R, Tuple(vars))
end

# handle generic kwargs vars
function Projection(R::AbstractLazyReachSet; vars)
    return Projection(R, Tuple(vars))
end

# ================================================================
# Time-shifted reach-set
# ================================================================

struct ShiftedReachSet{N, RT<:AbstractLazyReachSet{N}} <: AbstractLazyReachSet{N}
    R::RT
    t0::N
end

# getter functions
@inline time_shift(srs::ShiftedReachSet) = srs.t0

# time domain interface
@inline tstart(srs::ShiftedReachSet) = tstart(srs.R) + time_shift(srs)
@inline tend(srs::ShiftedReachSet) = tend(srs.R) + time_shift(srs)
@inline tspan(srs::ShiftedReachSet) = TimeInterval(tstart(srs), tend(srs))

# ================================================================
# Taylor model reach set
# ================================================================

"""
    AbstractTaylorModelReachSet{N}

Abstract type for all reach sets types that represent a Taylor model.

### Notes

The parameter `N` refers to the numerical type of the representation.
"""
abstract type AbstractTaylorModelReachSet{N} <: AbstractReachSet{N} end

using TaylorModels: TaylorModel1, TaylorN

"""
    TaylorModelReachSet{N} <: AbstractTaylorModelReachSet{N}

Taylor model reach-set represented as a vector taylor models in one variable
(namely, the "time" variable) whose coefficients are multivariate polynomials
(namely in the "space" variables).

### Notes

The parameter `N` refers to the numerical type of the representation.

In `TMJets`, the space variables are normalized to the interval `[-1, 1]`.
"""
struct TaylorModelReachSet{N} <: AbstractTaylorModelReachSet{N}
    X::Vector{TaylorModel1{TaylorN{N}, N}}
    Δt::TimeInterval
end

# interface functions
@inline set(R::TaylorModelReachSet) = R.X
@inline setrep(::TaylorModelReachSet{N}) where {N} = Vector{TaylorModel1{TaylorN{N}, N}}
@inline setrep(::Type{TaylorModelReachSet{N}}) where {N} = Vector{TaylorModel1{TaylorN{N}, N}}
@inline tstart(R::TaylorModelReachSet) = inf(R.Δt)
@inline tend(R::TaylorModelReachSet) = sup(R.Δt)
@inline tspan(R::TaylorModelReachSet) = R.Δt
@inline dim(R::TaylorModelReachSet) = get_numvars()
@inline vars(R::TaylorModelReachSet) = Tuple(Base.OneTo(length(R.X)),)

# useful constants
@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline unitBox(m) = IntervalBox(IA.Interval(0.0, 1.0), m)
@inline symBox(n::Integer) = IntervalBox(symI, n)
const zeroI = IA.Interval(0.0) # TODO use number type
const oneI = IA.Interval(1.0)
const symI = IA.Interval(-1.0, 1.0)

function shift(R::TaylorModelReachSet, t0::Number)
    return TaylorModelReachSet(set(R), tspan(R) + t0)
end

function reconstruct(R::TaylorModelReachSet, X)
    return TaylorModelReachSet(X, tspan(R))
end

function project(R::TaylorModelReachSet, vars::NTuple{D, M}) where {D, M<:Integer}
    throw(ArgumentError("the concrete projection of Taylor model reach-set is not " *
            "available; try first to overapproximate the Taylor model and the  project"))
end

# overapproximate taylor model reachset with one hyperrectangle
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle}) where {N}
    # dimension of the reachset
    D = dim(R)

    # pick the time domain of the given TM (same in all dimensions)
    t0 = tstart(R)
    Δt = tspan(R)
    # tdom defined below is the same as Δt - t0, but the domain inclusion check
    # in TM.evauate may fail, so we unpack the domain again here; also note that
    # by construction the TMs in time are centered at zero
    tdom = TM.domain(first(set(R)))

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = TM.evaluate(set(R), tdom)

    # evaluate the spatial variables in the symmetric box
    Bn = symBox(D)
    X̂ib = IntervalBox([TM.evaluate(X_Δt[i], Bn) for i in 1:D]...)
    X̂ = convert(Hyperrectangle, X̂ib)

    return ReachSet(X̂, Δt)
end

function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle}, nparts) where {N}
    # dimension of the reachset
    D = dim(R)

    # pick the time domain of the given TM (same in all dimensions)
    t0 = tstart(R)
    Δt = tspan(R)
    # tdom defined below is the same as Δt - t0, but the domain inclusion check
    # in TM.evauate may fail, so we unpack the domain again here; also note that
    # by construction the TMs in time are centered at zero
    tdom = TM.domain(first(set(R)))

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = TM.evaluate(set(R), tdom)

    # evaluate the spatial variables in the symmetric box
    partition = IA.mince(symBox(D), nparts)
    X̂ = Vector{Hyperrectangle{N, SVector{D, N}, SVector{D, N}}}(undef, length(partition))
    @inbounds for (i, Bi) in enumerate(partition)
        X̂ib = IntervalBox([TM.evaluate(X_Δt[i], Bi) for i in 1:D])
        X̂[i] = convert(Hyperrectangle, X̂ib)
    end
    #return ReachSet(UnionSetArray(X̂), Δt) # but UnionSetArray is not yet a lazyset
    return ReachSet(ConvexHullArray(X̂), Δt)
end

function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope}) where {N}
    # dimension of the reachset
    n = dim(R)

    # pick the time domain of the given TM (same in all dimensions)
    t0 = tstart(R)
    Δt = tspan(R)
    # tdom defined below is the same as Δt - t0, but the domain inclusion check
    # in TM.evauate may fail, so we unpack the domain again here; also note that
    # by construction the TMs in time are centered at zero
    tdom = TM.domain(first(set(R)))

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X = set(R)
    X_Δt = TM.evaluate(X, tdom)

    # builds the associated taylor model for each coordinate j = 1...n
    #  X̂ is a TaylorModelN whose coefficients are intervals
    X̂ = [TaylorModelN(X_Δt[j], X[j].rem, zeroBox(n), symBox(n)) for j in 1:n]

    # compute floating point rigorous polynomial approximation
    # fX̂ is a TaylorModelN whose coefficients are floats
    fX̂ = TaylorModels.fp_rpa.(X̂)

    # LazySets can overapproximate a Taylor model with a Zonotope
    Zi = overapproximate(fX̂, Zonotope)
    return ReachSet(Zi, Δt)
end

# ================================================================
# Template reach set
# ================================================================

using LazySets.Approximations: AbstractDirections

"""
    TemplateReachSet{N, VN, TN<:AbstractDirections{N, VN}, SN<:AbstractVector{N}} <: AbstractLazyReachSet{N}

Reach set that stores the support function of a set at a give set of directions.

### Notes

The struct has the following parameters:

- `N`  -- refers to the numerical type of the representation.
- `VN` -- refers to the vector type of the template
- `TN` -- refers to the template type
- `SN` -- vector type that holds the support function evaluations

Concrete subtypes of `AbstractDirections` are defined in the `LazySets` library.

This reach-set implicitly represents a set by a set of directions and support
functions. `set(R::TemplateReachSet)` returns a polyhedron in half-space
representation.

Apart from the getter functions inherited from the `AbstractReachSet` interface,
the following methods are available:

- `directions(R)`  -- return the set of directions normal to the faces of this reach-set
- `sup_func(R)`    -- return the vector of support function evaluations
- `sup_func(R, i)` -- return the `i`-th coordinate of the vector of support function evaluatons
"""
struct TemplateReachSet{N, VN, TN<:AbstractDirections{N, VN}, SN<:AbstractVector{N}} <: AbstractLazyReachSet{N}
    dirs::TN
    sf::SN # TODO fix sf::Vector{N} ?
    Δt::TimeInterval
end

# constructor from a given set (may overapproximate)
function TemplateReachSet(dirs::AbstractDirections, X::LazySet, Δt::TimeInterval)
    sfunc = [ρ(d, X) for d in dirs]
    return TemplateReachSet(dirs, sfunc, Δt)
end

# constructor from a given reach-set (may overapproximate)
function TemplateReachSet(dirs::AbstractDirections, R::AbstractLazyReachSet)
    return TemplateReachSet(dirs, set(R), tspan(R))
end

# implement abstract reachset interface
# TODO: use HPolyhedron or HPolytope if the set is bounded or not
set(R::TemplateReachSet) = HPolytope([HalfSpace(di, R.sf[i]) for (i, di) in enumerate(R.dirs)])
setrep(::Type{<:TemplateReachSet{N, VN}}) where {N, VN} = HPolyhedron{N, VN}
setrep(::TemplateReachSet{N, VN}) where {N, VN} = HPolyhedron{N, VN}
tspan(R::TemplateReachSet) = R.Δt
tstart(R::TemplateReachSet) = inf(R.Δt)
tend(R::TemplateReachSet) = sup(R.Δt)
dim(R::TemplateReachSet) = dim(R.dirs)
vars(R::TemplateReachSet) = Tuple(Base.OneTo(dim(R)),) # TODO rename to vars ?

directions(R::TemplateReachSet) = R.dirs # TODO rename to dirs ?
sup_func(R::TemplateReachSet) = R.sf # TODO rename?
sup_func(R::TemplateReachSet, i) = R.sf[i]

# overapproximate a given reach set with the template, concretely
#= TODO: remove?
function overapproximate(R::AbstractLazyReachSet, dirs::Vector{VN}) where {VN}
    Δt = tspan(R)
    sup_func = map(d -> ρ(d, R), dirs)
    return TemplateReachSet(dirs, sup_func, Δt)
end
=#

# convenience functions to get support directions of a given set
_getdirs(X::LazySet) = [c.a for c in constraints_list(X)]
LazySets.CustomDirections(X::LazySet) = CustomDirections(_getdirs(X), dim(X)) # TODO may use isboundedtype trait
