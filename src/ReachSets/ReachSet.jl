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
struct ReachSet{N,ST<:LazySet{N}} <: AbstractLazyReachSet{N}
    X::ST
    Δt::TimeInterval
end

# abstract reach set interface functions
set(R::ReachSet)                                          = R.X
setrep(R::ReachSet{N,ST}) where {N,ST<:LazySet{N}}        = ST
setrep(::Type{ReachSet{N,ST}}) where {N,ST<:LazySet{N}}   = ST
setrep(::AbstractVector{RT}) where {RT<:AbstractReachSet} = setrep(RT)
tstart(R::ReachSet)                                       = tstart(R.Δt)
tend(R::ReachSet)                                         = tend(R.Δt)
tspan(R::ReachSet)                                        = R.Δt
dim(R::ReachSet)                                          = dim(R.X)
vars(R::ReachSet)                                         = Tuple(Base.OneTo(dim(R.X)))

# no-op
function Base.convert(T::Type{ReachSet{N,ST}}, R::ReachSet{N,ST}) where {N,ST<:LazySet{N}}
    return R
end

# convert to a reach-set with another underlying set representation
function Base.convert(T::Type{ReachSet{N,LT}},
                      R::ReachSet{N,ST}) where {N,ST<:LazySet{N},LT<:LazySet{N}}
    X = convert(LT, set(R))
    return ReachSet(X, tspan(R))
end

function shift(R::ReachSet, t0::Number)
    return ReachSet(set(R), tspan(R) + t0)
end

function reconstruct(R::ReachSet, Y::LazySet)
    return ReachSet(Y, tspan(R))
end

# constructor with a time point
function ReachSet(X::ST, t::Real) where {N,ST<:LazySet{N}}
    return ReachSet(X, interval(t))
end

#=
# specialized concrete projection of zonotopic set by removing zero generators
function linear_map(R::ReachSet{N, <:AbstractZonotope{N}}, variables::NTuple{D,Int}; check_vars::Bool=true) where {D}
    πR = linear(R, variables, check_vars)
    X = remove_zero_generators(set(πR))
    SparseReachSet(proj, tspan(R), variables)
end

function project(R::ReachSet{N, <:AbstractZonotope{N}}, vars) where {N}
    πR = linear_map(M, R)
    dt = tspan(R)
    if isnothing(vars)
        X = remove_zero_generators(set(πR))
        v = ReachabilityAnalysis.vars(R)
    else
        Y = project(πR, vars)
        X = remove_zero_generators(set(Y))
        v = ReachabilityAnalysis.vars(Y)
    end
    return SparseReachSet(X, dt, v)
end
=#
