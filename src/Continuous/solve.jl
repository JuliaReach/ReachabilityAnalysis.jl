export solve

# ======================================
# Solve function for continuous systems
# ======================================

"""
    solve(ivp::IVP{<:AbstractContinuousSystem}, tspan, alg; kwargs...)

Solves the initial-value problem defined by `ivp` over the time span `tspan`,
using the algorihm `alg`. If no algorithm is given, a default algorithm is chosen.

### Input

- `ivp`   -- initial-value problem
- `tspan` -- time span for this initial-value problem
- `alg`   -- reachability algorithm

Additional options are passed as arguments or keyword arguments.
See the online documentation for examples.

### Output

A solution type.
"""
function solve(ivp::IVP{<:AbstractContinuousSystem}, args...; kwargs...)
    # preliminary checks
    _check_dim(ivp)

    # retrieve time span and continuous post operator algorithm
    tspan = _get_tspan(args...; kwargs...)
    cpost = _get_cpost(ivp, tspan, args...; kwargs...)

    # run the continuous-post operator
    F = post(cpost, ivp, tspan, args...; kwargs...)

    if haskey(kwargs, :save_traces) && (kwargs[:save_traces] == true)
        @requires DifferentialEquations
        error("saving traces is not implemented yet")
        # compute trajectories, cf. ensemble simulation
        # traces =
        # sol = ReachSolution(F, cpost, traces) # new solution type?
    else
        # wrap the flowpipe and algorithm in a solution structure
        sol = ReachSolution(F, cpost)
    end

    return sol
end

#=
function solve_continuous(ivp, cpost, tspan, args...; kwargs...)


    sol = post(ivp, cpost, tspan, args...; kwargs...)

    return sol
end
=#

# ==================
# Argument handling
# ==================

# extend dimension for common IntervalArithmetic types
LazySets.dim(::IA.Interval) = 1
LazySets.dim(::IA.IntervalBox{D, N}) where {D, N} = D

function _check_dim(ivp; throw_error::Bool=true)
    S = system(ivp)
    n = statedim(S)
    X0 = initial_state(ivp)
    if X0 isa LazySet
        d = dim(X0)
    elseif X0 isa AbstractVector
        d = length(X0)
    elseif X0 isa Number
        d = 1
    else
        throw(ArgumentError("the type of the initial condition, $(typeof(X0)), cannot be handled"))
    end

    if n == d
        return true
    else
        throw_error || throw(ArgumentError("the state-space dimension should match the " *
                                "dimension of the initial state, but they are of size " *
                                "$n and $(dim(X0)) respectively"))
        return false
    end
end

@inline _promote_tspan((t1, t2)::Tuple{T, T}) where {T} = (t1, t2)
@inline _promote_tspan((t1, t2)::Tuple{T, S}) where {T, S} = promote(t1, t2)
@inline _promote_tspan(tspan::Number) = (zero(tspan), tspan)
@inline _promote_tspan(tspan::IA.Interval) = (inf(tspan), sup(tspan))
@inline _promote_tspan(tspan::Interval) = (min(tspan), max(tspan))
@inline function _promote_tspan(tspan::AbstractArray)
    if length(tspan) == 2
        return (first(tspan), last(tspan))
    else
        throw(ArgumentError("the length of tspan must be two (and preferably, " *
                            "`tspan` should be a tuple, i.e. (0.0, 1.0)), but " *
                            "it is of length $(length(tspan))"))
    end
end

function _get_tspan(args...; kwargs...)
    got_tspan = haskey(kwargs, :tspan)
    got_T = haskey(kwargs, :T)

    if got_tspan && got_T
        throw(ArgumentError("cannot parse the time horizon `T` and the " *
            "time span `tspan` simultaneously; use one or the other"))
    end

    if got_tspan
        tspan = kwargs[:tspan]
        tspan = _promote_tspan(tspan)
    elseif got_T
        T = kwargs[:T]
        tspan = (zero(T), T)
    else
        throw(ArgumentError("the time span has not been specified, " *
            "but is required for `solve`; you should specify either the time horizon " *
            "`T=...` or the time span `tspan=...`"))
    end
    return TimeInterval(tspan[1], tspan[2])
end

# return the time horizon given a time span
# the check_positive flag is used for algorithms that do not support negative
# times
function _get_T(tspan::TimeInterval, check_zero::Bool=true, check_positive::Bool=true)
    t0 = inf(tspan)
    if check_zero
        @assert iszero(t0) "this algorithm can only handle zero initial time"
    end
    T = sup(tspan)
    if check_positive
        @assert T > 0 "the time horizon should be positive"
    end
    return T
end

function _get_cpost(ivp, tspan, args...; kwargs...)
    got_alg = haskey(kwargs, :alg)
    got_algorithm = haskey(kwargs, :algorithm)
    got_opC = haskey(kwargs, :opC)
    no_args = isempty(args) || args[1] === nothing

    if got_alg && no_args
        opC = kwargs[:alg]
    elseif got_algorithm && no_args
        opC = kwargs[:algorithm]
    elseif got_opC && no_args
        opC = kwargs[:opC]
    elseif !no_args && typeof(args[1]) <: AbstractContinuousPost
        opC = args[1]
    else
        opC = _default_cpost(ivp, tspan, args... ; kwargs...)
    end
    return opC
end

# ===================
# Algorithm defaults
# ===================

const DEFAULT_NSTEPS = 100

# If the system is affine, algorithm `GLGM06` is used.
# Otherwise, algorithm `TMJets` is used.
function _default_cpost(ivp::IVP{<:AbstractContinuousSystem}, tspan, args...; kwargs...)
    if isaffine(ivp)
        if haskey(kwargs, :δ)
            δ = kwargs[:δ]
        elseif haskey(kwargs, :N)
            N = kwargs[:N]
            δ = diam(tspan) / N
        else
            δ = diam(tspan) / DEFAULT_NSTEPS
        end
        opC = GLGM06(δ=δ)
    else
        opC = TMJets()
    end
    return opC
end
