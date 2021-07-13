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

# constructor with a time point
function SparseReachSet(X::ST, t::Real, vars::AbstractVector) where {N, ST<:LazySet{N}}
    return SparseReachSet(X, interval(t), vars)
end
