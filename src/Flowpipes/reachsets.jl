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
LazySets.dim(R::AbstractLazyReachSet) = dim(set(R))

LazySets.ρ(d::AbstractVector, R::AbstractLazyReachSet) = ρ(d, set(R))
LazySets.ρ(d::AbstractVector, R::Vector{<:AbstractLazyReachSet}) = map(Ri -> ρ(d, set(Ri)), R) |> maximum
LazySets.ρ(d::AbstractVector, R::SubArray{<:AbstractLazyReachSet}) = map(Ri -> ρ(d, set(Ri)), R) |> maximum

LazySets.σ(d::AbstractVector, R::AbstractLazyReachSet) = σ(d, set(R))
LazySets.σ(d::AbstractVector, R::Vector{<:AbstractLazyReachSet}) = _σ_vec(d, Rvec)
LazySets.σ(d::AbstractVector, R::SubArray{<:AbstractLazyReachSet}) = _σ_vec(d, Rvec)

function _σ_vec(d, Rvec)
    σarray = map(Ri -> σ(d, set(Ri)), Rvec)
    ρarray = map(vi -> dot(d, vi), σarray)
    m = argmax(ρarray)
    return σarray[m]
end

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

LazySets.isdisjoint(R::AbstractLazyReachSet, Y::LazySet) = isdisjoint(set(R), Y)

function LinearMap(M::Union{AbstractMatrix, Number}, R::AbstractLazyReachSet)
    return reconstruct(R, LinearMap(M, set(R)))
end

function linear_map(M::AbstractMatrix, R::AbstractLazyReachSet)
    return reconstruct(R, linear_map(M, set(R)))
end

function linear_map(M::AbstractMatrix, R::Vector{<:AbstractLazyReachSet})
    return map(Ri -> reconstruct(Ri, linear_map(M, set(Ri))), R)
end

function overapproximate(R::AbstractLazyReachSet, func)
    return reconstruct(R, overapproximate(set(R), func))
end

# concrete projection extensions
project(R::AbstractLazyReachSet, vars::AbstractVector{M}) where {M<:Integer} = project(R, Tuple(vars))
project(R::AbstractLazyReachSet; vars) = project(R, Tuple(vars))

# projectin of an array of reach-sets
_project_vec(R, vars) = map(Ri -> project(Ri, Tuple(vars)), R)
project(R::Vector{<:AbstractLazyReachSet}, vars::AbstractVector{M}) where {M<:Integer} = _project_vec(R, vars)
project(R::Vector{<:AbstractLazyReachSet}; vars) = _project_vec(R, vars)
project(R::SubArray{<:AbstractLazyReachSet}, vars) = _project_vec(R, vars)
project(R::SubArray{<:AbstractLazyReachSet}; vars) = _project_vec(R, vars)

function project(R::AbstractLazyReachSet, M::AbstractMatrix; vars=nothing)
    πR = linear_map(M, R)
    return isnothing(vars) ? πR : project(πR, vars)
end

function project(R::Vector{<:AbstractLazyReachSet}, M::AbstractMatrix; vars=nothing)
    πR = map(Ri -> linear_map(M, Ri), R)
    return isnothing(vars) ? πR : project(πR, vars)
end

# membership test
function ∈(x::AbstractVector{N}, R::AbstractLazyReachSet{N}) where {N}
    return ∈(x, set(R))
end

# splitting a reach-set according to a given partition; the partition should be
# a vector of integers
function LazySets.split(R::AbstractLazyReachSet, partition)
    Y = split(set(R), partition)
    [reconstruct(R, y) for y in Y]
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
    project(R::AbstractLazyReachSet, variables::NTuple{D, M};
            check_vars::Bool=true) where {D, M<:Integer}

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

    vR = vars(R)
    vRvec = collect(vR)

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(variables, 0) ⊆ vR)
        throw(ArgumentError("the variables $vars do not belong to the variables " *
                            " of this reach-set, $(vR)"))
    end

    if 0 ∈ variables  # the projection involves "time"
        vars_idx = _get_vars_idx(variables, vcat(0, vRvec))
        Δt = convert(Interval, tspan(R))
        proj =  _project(Δt × set(R), vars_idx)
    else
        vars_idx = _get_vars_idx(variables, vRvec)
        proj = _project(set(R), vars_idx)
    end

    return SparseReachSet(proj, tspan(R), variables)
end

# assumes that variables is a subset of all_variables
@inline function _get_vars_idx(variables, all_variables)
    vars_idx = indexin(variables, all_variables) |> Vector{Int}
end

# concrete projection onto a single variable
function project(R::AbstractLazyReachSet, variable::Int; check_vars::Bool=true)
    return project(R, (variable,), check_vars=check_vars)
end

# lazy projection of a reach-set
function Projection(R::AbstractLazyReachSet, variables::NTuple{D, M},
                    check_vars::Bool=true) where {D, M<:Integer}

    vR = vars(R)
    vRvec = collect(vR)

    # TODO: make vars check faster, specific for ReachSets and number of vars D
    if check_vars && !(setdiff(variables, 0) ⊆ vR)
        throw(ArgumentError("the variables $variables do not belong to the variables " *
                            " of this reach-set, $(vR)"))
    end

    if 0 ∈ variables  # the projection involves "time"
        vars_idx = _get_vars_idx(variables, vcat(0, vRvec))
        Δt = convert(Interval, tspan(R))
        proj =  _Projection(Δt × set(R), vars_idx)
    else
        vars_idx = _get_vars_idx(variables, vRvec)
        proj = _Projection(set(R), vars_idx)
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

Taylor model reach-set represented as a vector of taylor models in one variable
(namely, the "time" variable) whose coefficients are multivariate polynomials
(namely in the "space" variables). It is assumed that the time domain is the same
for all components.

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

# overload getter functions for the taylor model
# we assume that the first element is representative
domain(R::TaylorModelReachSet) = domain(first(R.X))
remainder(R::TaylorModelReachSet) = remainder.(R.X)
polynomial(R::TaylorModelReachSet) = polynomial.(R.X)
get_order(R::TaylorModelReachSet) = get_order.(R.X)
expansion_point(R::TaylorModelReachSet) = [Xi.x0 for Xi in R.X]

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

- `directions(R)`           -- return the set of directions normal to the faces of this reach-set
- `support_functions(R)`    -- return the vector of support function evaluations
- `support_functions(R, i)` -- return the `i`-th coordinate of the vector of support function evaluatons
"""
struct TemplateReachSet{N, VN, TN<:AbstractDirections{N, VN}, SN<:AbstractVector{N}} <: AbstractLazyReachSet{N}
    dirs::TN
    sf::SN
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

# abstract reachset interface
function set(R::TemplateReachSet)
    T = isbounding(R.dirs) ? HPolytope : HPolyhedron
    return T([HalfSpace(di, R.sf[i]) for (i, di) in enumerate(R.dirs)])
end

# FIXME requires adding boundedness property as a type parameter
setrep(::Type{<:TemplateReachSet{N, VN}}) where {N, VN} = HPolyhedron{N, VN}
setrep(::TemplateReachSet{N, VN}) where {N, VN} = HPolyhedron{N, VN}

tspan(R::TemplateReachSet) = R.Δt
tstart(R::TemplateReachSet) = inf(R.Δt)
tend(R::TemplateReachSet) = sup(R.Δt)
dim(R::TemplateReachSet) = dim(R.dirs)
vars(R::TemplateReachSet) = Tuple(Base.OneTo(dim(R)),)

directions(R::TemplateReachSet) = R.dirs
support_functions(R::TemplateReachSet) = R.sf
support_functions(R::TemplateReachSet, i::Int) = R.sf[i]

_collect(dirs::AbstractDirections) = collect(dirs)
_collect(dirs::CustomDirections) = dirs.directions

# overapproximate a given reach-set with a template
function overapproximate(R::AbstractLazyReachSet, dirs::AbstractDirections) where {VN}
    return TemplateReachSet(dirs, R)
end
function overapproximate(R::AbstractLazyReachSet, dirs::Vector{VN}) where {VN}
    return TemplateReachSet(CustomDirections(dirs), R)
end

# concrete projection of a TemplateReachSet for a given direction
function project(R::TemplateReachSet, dir::AbstractVector{<:AbstractFloat}; vars=nothing)
    ρdir = ρdir₋ = nothing
    dir₋ = -dir
    @inbounds for (i, d) in enumerate(R.dirs)
        if iszero(norm(dir - d))
            ρdir = R.sf[i]
        elseif iszero(norm(dir₋ - d))
            ρdir₋ = -R.sf[i]
        end
    end
    if isnothing(ρdir)
        ρdir = ρ(dir, R)
    end
    if isnothing(ρdir₋)
        ρdir₋ = -ρ(ρdir₋, R)
    end
    πR = ReachSet(Interval(min(ρdir, ρdir₋), max(ρdir, ρdir₋)), tspan(R))
    return isnothing(vars) ? πR : project(πR, vars)
end

# concrete projection of a vector of TemplateReachSet for a given direction
function project(Xk::Vector{<:TemplateReachSet}, dir::AbstractVector{<:AbstractFloat}; vars=nothing)
    return map(X -> project(X, dir, vars=vars), Xk)
end

# convenience functions to get support directions of a given set
_getdirs(X::LazySet) = [c.a for c in constraints_list(X)]
LazySets.CustomDirections(X::LazySet) = CustomDirections(_getdirs(X), dim(X)) # TODO may use isboundedtype trait
