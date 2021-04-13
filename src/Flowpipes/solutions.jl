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

TODO

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

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`AbstractReachSet`](@ref)s
- `alg`      -- algorihm used
- `options`  -- the dictionary of options
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
flowpipe(sol::ReachSolution) = flowpipe(sol.F)
tstart(sol::ReachSolution) = tstart(sol.F)
tend(sol::ReachSolution) = tend(sol.F)
tspan(sol::ReachSolution) = tspan(sol.F)
tstart(sol::ReachSolution, arr::AbstractVector) = tstart(sol.F, arr)
tend(sol::ReachSolution, arr::AbstractVector) = tend(sol.F, arr)
tspan(sol::ReachSolution, arr::AbstractVector) = tspan(sol.F, arr)
ensemble(sol::ReachSolution) = get(sol.ext, :ensemble, nothing)

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

# evaluation interface
Base.getindex(sol::ReachSolution, t::Float64) = getindex(sol.F, t)
(sol::ReachSolution)(t::Float64) = sol.F(t)
(sol::ReachSolution)(t::Number) = sol.F(t)
(sol::ReachSolution)(dt::IA.Interval{Float64}) = sol.F(dt)

function overapproximate(sol::ReachSolution{FT}, args...) where {FT<:AbstractFlowpipe}
    return ReachSolution(overapproximate(sol.F, args...), sol.alg, sol.ext)
end

function project(sol::ReachSolution{FT}, args...) where {FT<:AbstractFlowpipe}
    return project(sol.F, args...)
end

# convenience alias to match the usage in the plot recipe
function project(sol::ReachSolution{FT}; vars) where {FT<:AbstractFlowpipe}
    return project(sol.F, Tuple(vars))
end

function shift(sol::ReachSolution{<:AbstractFlowpipe}, t0::Number)
    return ReachSolution(shift(sol.F, t0), sol.alg, sol.ext)
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

# -------------------------------------------
# Specialized methods for template solutions
# -------------------------------------------

function support_function_matrix(sol::ReachSolution{<:Flowpipe{N, <:TemplateReachSet}}) where {N}
    return support_function_matrix(flowpipe(sol))
end

function flatten(sol::ReachSolution{<:Flowpipe{N, <:TemplateReachSet}}, rows=(1, 2)) where {N}
    return flatten(flowpipe(sol), rows)
end
