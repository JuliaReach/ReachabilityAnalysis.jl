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

- `linear_map`  -- concrete linear map of a reach-set
- `project`     -- projection of a reach-set
- `shift`       -- time-shift of a reach-set
- `vars`        -- tuple of integers associated to the variables of the given reach-set
"""
abstract type AbstractLazyReachSet{N} <: AbstractReachSet{N} end

# Implement LazySets interface
LazySets.dim(R::AbstractLazyReachSet) = dim(set(R))

LazySets.ρ(d::AbstractVector, R::AbstractLazyReachSet) = ρ(d, set(R))
function LazySets.ρ(d::AbstractVector, R::Vector{<:AbstractLazyReachSet})
    return maximum(map(Ri -> ρ(d, set(Ri)), R))
end
function LazySets.ρ(d::AbstractVector, R::SubArray{<:AbstractLazyReachSet})
    return maximum(map(Ri -> ρ(d, set(Ri)), R))
end

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
for ST in (LazySet, UnionSet, UnionSetArray)
    @eval Base.:⊆(R::AbstractLazyReachSet, X::$(ST)) = ⊆(set(R), X)
    @eval Base.:⊆(X::$(ST), R::AbstractLazyReachSet) = ⊆(X, set(R))
end
Base.:⊆(R::AbstractLazyReachSet, S::AbstractLazyReachSet) = ⊆(set(R), set(S))
LazySets.area(R::AbstractLazyReachSet) = area(set(R))
LazySets.volume(R::AbstractLazyReachSet) = volume(set(R))
Base.convert(::Type{ST}, R::AbstractLazyReachSet) where {ST<:LazySet} = convert(ST, set(R))
Base.convert(::Type{<:IntervalBox}, R::AbstractLazyReachSet) = convert(IntervalBox, set(R))
LazySets.complement(R::AbstractLazyReachSet) = reconstruct(R, complement(set(R)))
Base.isempty(R::AbstractLazyReachSet) = isempty(set(R))

function LinearMap(M::Union{AbstractMatrix,Number}, R::AbstractLazyReachSet)
    return reconstruct(R, LinearMap(M, set(R)))
end

Base.:(*)(M::AbstractMatrix, R::AbstractLazyReachSet) = LinearMap(M, R)
Base.:(*)(α::AbstractFloat, R::AbstractLazyReachSet) = reconstruct(R, α * set(R))

function linear_map(M::AbstractMatrix, R::AbstractLazyReachSet)
    return reconstruct(R, linear_map(M, set(R)))
end

function linear_map(M::AbstractMatrix, R::Vector{<:AbstractLazyReachSet})
    return map(Ri -> reconstruct(Ri, linear_map(M, set(Ri))), R)
end

function overapproximate(R::AbstractLazyReachSet, func)
    return reconstruct(R, overapproximate(set(R), func))
end

# ------------------------------------
# Concrete projection of a reach-set
# ------------------------------------

"""
    project(R::AbstractLazyReachSet, variables::NTuple{D,Int};
            check_vars::Bool=true) where {D}

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
function project(R::AbstractLazyReachSet, variables::NTuple{D,Int};
                 check_vars::Bool=true) where {D}
    vR = vars(R)
    vRvec = collect(vR)

    if check_vars && !(setdiff(variables, 0) ⊆ vR)
        throw(ArgumentError("the variables $vars do not belong to the variables " *
                            " of this reach-set, $(vR)"))
    end

    if 0 ∈ variables  # the projection involves "time"
        vars_idx = _get_vars_idx(variables, vcat(0, vRvec))
        Δt = convert(Interval, tspan(R))
        proj = project(Δt × set(R), vars_idx)
    else
        vars_idx = _get_vars_idx(variables, vRvec)
        proj = project(set(R), vars_idx)
    end

    return SparseReachSet(proj, tspan(R), variables)
end

# assumes that variables is a subset of all_variables
@inline function _get_vars_idx(variables, all_variables)
    return vars_idx = Vector{Int}(indexin(variables, all_variables))
end

# concrete projection onto a single variable
function project(R::AbstractLazyReachSet, variable::Int; check_vars::Bool=true)
    return project(R, (variable,); check_vars=check_vars)
end

# concrete projection overloads
function project(R::AbstractLazyReachSet, vars::AbstractVector{Int})
    return project(R, Tuple(vars))
end
project(R::AbstractLazyReachSet; vars) = project(R, Tuple(vars))

# projection of an array of reach-sets
_project_vec(R, vars) = map(Ri -> project(Ri, Tuple(vars)), R)
project(R::AbstractVector{<:AbstractLazyReachSet}, vars::VecOrTupleOrInt) = _project_vec(R, vars)
project(R::AbstractVector{<:AbstractLazyReachSet}; vars::VecOrTupleOrInt) = _project_vec(R, vars)

# -------------------------------
# Lazy projection of a reach-set
# -------------------------------

function Projection(R::AbstractLazyReachSet, variables::NTuple{D,Int},
                    check_vars::Bool=true) where {D}
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
        proj = Projection(Δt × set(R), vars_idx)
    else
        vars_idx = _get_vars_idx(variables, vRvec)
        proj = Projection(set(R), vars_idx)
    end

    return SparseReachSet(proj, tspan(R), variables)
end

# handle generic vars vector
function Projection(R::AbstractLazyReachSet, vars::AbstractVector{Int})
    return Projection(R, Tuple(vars))
end

# handle generic kwargs vars
function Projection(R::AbstractLazyReachSet; vars)
    return Projection(R, Tuple(vars))
end

# ----------------------
# Additional extensions
# ----------------------

# membership test
function ∈(x::AbstractVector{N}, R::AbstractLazyReachSet{N}) where {N}
    return ∈(x, set(R))
end

# splitting a reach-set according to a given partition; the partition should be
# a vector of integers
function LazySets.split(R::AbstractLazyReachSet, partition)
    Y = split(set(R), partition)
    return [reconstruct(R, y) for y in Y]
end
