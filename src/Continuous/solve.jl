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

Additional options are passed as arguments or keyword arguments; see the notes
below for details. See the online documentation for examples.

### Output

The solution of a reachability problem, as an instance of a `ReachSolution`.

### Notes

- Use the `alg`, `algorithm` or `opC` keyword arguments to specify the algorithm
  to solve the initial-value problem. Algorithm-specific options should be passed
  to the algorithm constructor as well.

- Use the `tspan` keyword argument to specifying the time span; it can be:
    - a tuple,
    - an interval, or
    - a vector with two components.

- Use the `T` keyword argument to specify the time horizon; the initial time is
  then assumed to be zero.

- Use the `static` keyword argument to force conversion to static arrays in the
  algorithm (should be supported by the algorithm).

- Use the `NSTEPS` keyword argument to specify the number of discrete steps
  solved in the set-based recurrence.

- Use the `threading` option to use multi-threading parallelism. This option applies
  for initial-value problems whose initial condition is a vector of sets.
"""
function solve(ivp::IVP{<:AbstractContinuousSystem}, args...; kwargs...)
    # preliminary checks
    _check_dim(ivp)

    # get time span (or the emptyset if NSTEPS was specified)
    tspan = _get_tspan(args...; kwargs...)

    # get the continuous post or find a default one
    cpost = _get_cpost(ivp, args...; kwargs...)
    if cpost == nothing
        cpost = _default_cpost(ivp, tspan; kwargs...)
    end

    # run the continuous-post operator
    F = post(cpost, ivp, tspan; kwargs...)

    if haskey(kwargs, :save_traces) && kwargs[:save_traces]
        @requires DifferentialEquations
        error("saving traces is not implemented yet")
        # compute trajectories, cf. ensemble simulation
        # traces = ...
        # sol = ReachSolution(F, cpost, traces) # new solution type?
    else
        # wrap the flowpipe and algorithm in a solution structure
        sol = ReachSolution(F, cpost)
    end

    return sol
end

# solve for distributed initial conditions; uses multi-threaded implementation by default
function solve(ivp::IVP{AT, VT}, args...; kwargs...) where {AT<:AbstractContinuousSystem,
                                                            ST<:LazySet, VT<:AbstractVector{ST}}
    _check_dim(ivp)
    tspan = _get_tspan(args...; kwargs...)
    cpost = _get_cpost(ivp, args...; kwargs...)
    if cpost == nothing
        cpost = _default_cpost(ivp, tspan; kwargs...)
    end

    X0 = initial_state(ivp)
    S = system(ivp)
    # accept "threading"
    threading = haskey(kwargs, :threading) ? kwargs[:threading] : false # TEMP (Threads.nthreads() > 1)

    F = _solve_distributed(cpost, S, X0, tspan, Val(threading); kwargs...)
    return ReachSolution(MixedFlowpipe(F), cpost)
end

function _solve_distributed(cpost, S, X0, tspan, threading::Val{false}; kwargs...)
    return [post(cpost, IVP(S, X0i), tspan; kwargs...) for X0i in X0]
end

function _solve_distributed(cpost, S, X0, tspan, threading::Val{true}; kwargs...)
    nsets = length(X0)
    FT = Flowpipe{numtype(cpost), rsetrep(cpost)}
    sol_tot = Vector{FT}(undef, nsets)

    Threads.@threads for i in 1:length(X0)
        sol_i = ReachabilityAnalysis.post(cpost, IVP(S, X0[i]), tspan; kwargs...)
        sol_tot[i] = sol_i
    end
    return sol_tot
end

# ====================================
# Continuous post operator interface
# ====================================

"""
    AbstractPost

Abstract supertype of all post operator types.
"""
abstract type AbstractPost end

"""
    AbstractContinuousPost

Abstract supertype of all continuous post operators.
"""
abstract type AbstractContinuousPost <: AbstractPost end

setrep(::InitialValueProblem{HS, ST}) where {HS, ST} = ST

# ==================
# Argument handling
# ==================

# extend dimension for common IntervalArithmetic types
LazySets.dim(::IA.Interval) = 1
LazySets.dim(::IA.IntervalBox{D, N}) where {D, N} = D

# lazy sets (or sets that behave like such -- TODO update once UnionSet <: LazySet)
_dim(X::Union{<:LazySet, <:IA.Interval, <:IA.IntervalBox, <:UnionSet, <:UnionSetArray}) = dim(X)

# singleton elements
_dim(X::Number) = 1
_dim(X::AbstractVector{N}) where {N<:Number} = length(X)

# vector of sets
function  _dim(X::AbstractVector{UT}) where {UT<:Union{<:LazySet,
                <:IA.Interval, <:IA.IntervalBox, <:UnionSet, <:UnionSetArray}}
    n = _dim(first(X))
    all(X -> _dim(X) == n, X)  || throw(ArgumentError("dimension mismatch between " *
            "the initial sets in this array; expected only sets of dimension $n"))
    return n
end

_dim(X) = throw(ArgumentError("the type of the initial condition, $(typeof(X)), cannot be handled"))

function _check_dim(S, X0; throw_error::Bool=true)
    n = statedim(S)
    d = _dim(X0)
    n == d && return true

    if throw_error
        throw(ArgumentError("the state-space dimension should match " *
                            "the dimension of the initial state, but they are " *
                            "$n and $(dim(X0)) respectively"))
    end
    return false
end

function _check_dim(ivp::InitialValueProblem{<:AbstractContinuousSystem}; throw_error::Bool=true)
    S = system(ivp)
    X0 = initial_state(ivp)
    _check_dim(S, X0, throw_error=throw_error)
end

# promotion from tuple-like arguments
@inline _promote_tspan((t1, t2)::Tuple{T, T}) where {T} = TimeInterval(t1, t2)
@inline _promote_tspan((t1, t2)::Tuple{T, S}) where {T, S} = TimeInterval(promote(t1, t2))

# no-op, corresponds to (inf(tspan), sup(tspan))
@inline _promote_tspan(tspan::IA.Interval) = tspan

# no-op, takes interval wrapped data; corresponds to (min(tspan), max(tspan))
@inline _promote_tspan(tspan::Interval) = tspan.dat

# number T defaults to a time interval of the form [0, T]
@inline _promote_tspan(tspan::Number) = TimeInterval(zero(tspan), tspan)

# catch-all array like tspan
@inline function _promote_tspan(tspan::AbstractArray)
    @assert length(tspan) == 2 throw(ArgumentError("the length of `tspan` must be two, but " *
                                                   "it is of length $(length(tspan))"))
    return TimeInterval(tspan[1], tspan[2])
end

function _get_tspan(args...; kwargs...)
    got_tspan = haskey(kwargs, :tspan)
    got_T = haskey(kwargs, :T)
    got_NSTEPS = haskey(kwargs, :NSTEPS)

    # throw error for repeated arguments
    if got_tspan && got_T
        throw(ArgumentError("cannot parse the time horizon `T` and the " *
            "time span `tspan` simultaneously; use one or the other"))
    elseif got_tspan && got_NSTEPS
        throw(ArgumentError("cannot parse the time span `tspan` and the " *
            "number of steps `NSTEPS` simultaneously; use one or the other"))
    elseif got_T && got_NSTEPS
        throw(ArgumentError("cannot parse the time horizon `T` and the " *
            "number of steps `NSTEPS` simultaneously; use one or the other"))
    end

    # parse time span
    if got_tspan
        tspan =_promote_tspan(kwargs[:tspan])

    elseif got_T
        tspan = _promote_tspan(kwargs[:T])

    elseif got_NSTEPS
        tspan = emptyinterval() # defined a posteriori

    elseif !isempty(args) && applicable(_promote_tspan, args[1])
        # applies to tuples, vectors, LazySets.Interval, IA.Interval, and numbers
        tspan = _promote_tspan(args[1])
#=
        (args[1] isa Tuple{Float64, Float64} || # time span given as tuple
                              args[1] isa Vector{Float64}         || # time span given as vector
                              args[1] isa Interval ||  # time span given as LazySets.Interval
                              args[1] isa IntervalArithmetic.Interval || # time span given as IA.Interval
                              args[1] isa Real) # got time horizon as first argument
=#

    elseif  haskey(kwargs, :max_jumps)
        # got maximum number of jumps, applicable to hybrid systems
        tspan = emptyinterval()

    else

        # couldn't find time span => error
        throw(ArgumentError("the time span has not been specified, but is required " *
                            "for `solve`; you should specify either the time horizon " *
                            "`T=...`, the time span `tspan=...`, or the number of steps, `NSTEPS=...`"))
    end
    return tspan
end

# return the time horizon given a time span
# the check_positive flag is used for algorithms that do not support negative
# times
function _get_T(tspan::TimeInterval; check_zero::Bool=true, check_positive::Bool=true)
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

tstart(Δt::TimeInterval) = inf(Δt)
tend(Δt::TimeInterval) = sup(Δt)
tspan(Δt::TimeInterval) = Δt

function _get_cpost(ivp, args...; kwargs...)
    got_alg = haskey(kwargs, :alg)
    got_algorithm = haskey(kwargs, :algorithm)
    got_opC = haskey(kwargs, :opC)
    no_args = isempty(args) || args[1] === nothing

    # continous post was specified
    if got_alg
        cpost = kwargs[:alg]
    elseif got_algorithm
        cpost = kwargs[:algorithm]
    elseif got_opC
        cpost = kwargs[:opC]

    # check if either args[1] or args[2] are the post
    elseif !no_args && typeof(args[1]) <: AbstractContinuousPost
        cpost = args[1]
    elseif !no_args && length(args) == 2 && typeof(args[2]) <: AbstractContinuousPost
        cpost = args[2]
    elseif !no_args && length(args) > 2
        throw(ArgumentError("the number of arguments, $(length(args)), is not valid"))

    # no algorithm found => compute default
    else
        cpost = nothing
    end
    return cpost
end

# ===================
# Algorithm defaults
# ===================

# number of discrete steps used if only the time horizonis given
const DEFAULT_NSTEPS = 100

function _default_cpost(ivp::AbstractContinuousSystem, tspan; kwargs...)
    if isaffine(ivp)
        if haskey(kwargs, :δ)
            δ = kwargs[:δ]
        elseif haskey(kwargs, :N)
            N = kwargs[:N]
            δ = diam(tspan) / N
        elseif haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
            δ = diam(tspan) / NSTEPS
        elseif haskey(kwargs, :num_steps)
            num_steps = kwargs[:num_steps]
            δ = diam(tspan) / num_steps
        else
            δ = diam(tspan) / DEFAULT_NSTEPS
        end
        if statedim(ivp) == 1
            opC = INT(δ=δ)
        else
            static = haskey(kwargs, :static) ? kwargs[:static] : false
            # TODO pass order and dimension options as well
            opC = GLGM06(δ=δ, static=static)
        end
    else
        opC = TMJets()
    end
    return opC
end

# TODO refactor / update docstring

"""
    _default_cpost(ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)

### Input

- `ivp`   -- initial-value problem
- `tspan` -- time-span

### Output

A continuous post-operator that can handle the given initial-value problem.

### Notes

If the system is affine, then:

- If it is one-dimensional, algorithm `INT` is used, otherwise,
- Algorithm `GLGM06` is used.

If the system is not affine, then the algorithm `TMJets` is used.
"""
function _default_cpost(ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)
    _default_cpost(system(ivp), tspan; kwargs...)
end

# return a default algorithm, assuming that the first mode is representative
function _default_cpost(ivp::IVP{<:AbstractHybridSystem}, tspan; kwargs...)
    first_mode = mode(system(ivp), 1)
    _default_cpost(first_mode, tspan; kwargs...)
end
