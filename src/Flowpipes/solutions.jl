# ================================
# Abstract types
# ================================

"""
    AbstractSolution

Abstract supertype of all solution types of a rechability problem.
"""
abstract type AbstractSolution end

# ================================
# Property checking problem
# ================================

"""
    CheckSolution{T} <: AbstractSolution

Type that wraps the solution of a verification problem.

### Fields

property::PT
satisfied::Bool
vidx::Int
vtspan::TimeInterval
alg::ST

### Notes

This type contains the answer if the property is satisfied, and if not, it
contains the index at which the property might be violated for the first time.
"""
struct CheckSolution{PT, ST<:AbstractPost} <: AbstractSolution
    property::PT
    satisfied::Bool
    vidx::Int
    vtspan::TimeInterval
    alg::ST
end

property(sol::CheckSolution) = sol.property
satisfied(sol::CheckSolution) = sol.satisfied
violation_index(sol::CheckSolution) = sol.vidx
violation_tspan(sol::CheckSolution) = sol.vtspan
alg(sol::CheckSolution) = sol.alg

# ================================
# Reachability problem
# ================================

"""
    ReachSolution{FT<:AbstractFlowpipe, ST<:AbstractPost} <: AbstractSolution

Type that wraps the solution of a reachability problem as a flowpipe, the algorithm
used to obtain the flowpipe, and a dictionary to store additional data.

### Fields

- `F`    -- flowpipe
- `alg`  -- algorithm
- `ext`  -- extension dictionary to store additional data / options
"""
struct ReachSolution{FT<:AbstractFlowpipe, ST<:AbstractPost} <: AbstractSolution
    F::FT
    alg::ST
    ext::Dict{Symbol, Any} # dictionary used by extensions
end

# constructor from empty extension dictionary
function ReachSolution(F::FT, alg::ST) where {FT<:AbstractFlowpipe, ST<:AbstractPost}
    return ReachSolution(F, alg, Dict{Symbol, Any}())
end

# getter functions
flowpipe(sol::ReachSolution) = sol.F
tstart(sol::ReachSolution) = tstart(sol.F)
tend(sol::ReachSolution) = tend(sol.F)
tspan(sol::ReachSolution) = tspan(sol.F)
tstart(sol::ReachSolution, arr::AbstractVector) = tstart(sol.F, arr)
tend(sol::ReachSolution, arr::AbstractVector) = tend(sol.F, arr)
tspan(sol::ReachSolution, arr::AbstractVector) = tspan(sol.F, arr)

# NOTE: using sol.alg for setrep and rsetrep may be inaccurate as solutions
# do not change their sol.alg (eg. overapproximate(sol, Zonotope)) where sol is
# a Taylor model flowpipe keeps sol.alg being TMJets => we have to check
# for the rep's of the flowpipe
setrep(sol::ReachSolution{FT, ST}) where {FT, ST} = setrep(FT)
rsetrep(sol::ReachSolution{FT, ST}) where {FT, ST} = rsetrep(FT)
vars(sol::ReachSolution) = vars(sol.F)
numrsets(sol::ReachSolution) = numrsets(sol.F)

# iteration and indexing iterator interface
array(sol::ReachSolution) = array(sol.F)
Base.iterate(sol::ReachSolution) = iterate(sol.F)
Base.iterate(sol::ReachSolution, state) = iterate(sol.F, state)
Base.length(sol::ReachSolution) = length(sol.F)
Base.first(sol::ReachSolution) = first(sol.F)
Base.last(sol::ReachSolution) = last(sol.F)
Base.firstindex(sol::ReachSolution) = 1
Base.lastindex(sol::ReachSolution) = length(sol.F)
Base.getindex(sol::ReachSolution, i::Int) = getindex(sol.F, i)
Base.getindex(sol::ReachSolution, i::Number) = getindex(sol.F, i)
Base.getindex(sol::ReachSolution, I::AbstractVector) = getindex(sol.F, I)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(sol::ReachSolution, I::Int...) = getindex(sol.F, I...)
Base.eachindex(sol::ReachSolution) = eachindex(sol.F)
Base.keys(sol::ReachSolution) = keys(sol.F)

# other Base functions
Base.copy(sol::ReachSolution) = deepcopy(sol)

# evaluation interface
Base.getindex(sol::ReachSolution, t::Float64) = getindex(sol.F, t)
(sol::ReachSolution)(t::Float64) = sol.F(t)
(sol::ReachSolution)(t::Number) = sol.F(t)
(sol::ReachSolution)(dt::IA.Interval{Float64}) = sol.F(dt)

# common set interface
function overapproximate(sol::ReachSolution, args...)
    return overapproximate(sol.F, args...)
end

# concrete projection of a solution
project(sol::ReachSolution, args...) = project(sol.F, args...)

# convenience alias to match the usage in the plot recipe
project(sol::ReachSolution; vars) = project(sol.F, Tuple(vars))

# concrete projection given a projection matrix
function project(sol::ReachSolution, M::AbstractMatrix; vars=nothing)
    return project(sol.F, M; vars=vars)
end

# concrete projection of a solution for a given direction
function project(sol::ReachSolution, dir::AbstractVector{<:AbstractFloat}; vars=nothing)
    return project(sol.F, dir; vars=vars)
end

# lazy projection of a solution
Projection(sol::ReachSolution{<:AbstractFlowpipe}, vars) = Projection(sol.F, vars)
Projection(sol::ReachSolution{<:AbstractFlowpipe}; vars) = Projection(sol.F, vars)

function shift(sol::ReachSolution{<:AbstractFlowpipe}, t0::Number)
    return shift(sol.F, t0)
end

# LazySets interface falls back to the associated flowpipe
LazySets.dim(sol::ReachSolution) = dim(sol.F)
LazySets.ρ(d, sol::ReachSolution{FT}) where {FT<:AbstractFlowpipe} = ρ(d, sol.F)
LazySets.σ(d, sol::ReachSolution{FT}) where {FT<:AbstractFlowpipe} = σ(d, sol.F)

# further setops functions acting on solutions' flowpipes
∈(x::AbstractVector, sol::ReachSolution) = ∈(x, sol.F)

LazySets.is_intersection_empty(sol::ReachSolution, Y::LazySet) = is_intersection_empty(sol.F, Y)
LazySets.is_intersection_empty(sol::ReachSolution, Y::AbstractLazyReachSet) = is_intersection_empty(sol.F, set(Y))
linear_map(M::AbstractMatrix, sol::ReachSolution) = linear_map(M, sol.F)

# inclusion checks
Base.:⊆(sol::ReachabilityAnalysis.ReachSolution, X::LazySet) = ⊆(sol.F, X)
Base.:⊆(sol::ReachabilityAnalysis.ReachSolution, Y::AbstractLazyReachSet) = ⊆(sol.F, set(Y))
