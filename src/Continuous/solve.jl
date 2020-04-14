# ======================================
# Solve function for continuous systems
# ======================================

"""
    solve(ivp::IVP, tspan, alg; kwargs...)

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
"""
function solve(ivp::IVP, args...; kwargs...)
    # preliminary checks
    _check_dim(ivp)

    # retrieve time span and continuous post operator algorithm
    # NOTE: tspan is of the form (0, 0) if NSTEPS was specified
    tspan = _get_tspan(args...; kwargs...)
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

# solve for distributed initial conditions
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
    F = [post(cpost, IVP(S, X0i), tspan; kwargs...) for X0i in X0] |> MixedFlowpipe
    return ReachSolution(F, cpost)
end

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
    elseif X0 isa AbstractVector && eltype(X0) <: Number
        d = length(X0)
    elseif X0 isa AbstractVector && eltype(X0) <: LazySet
        d = dim(first(X0)) # assumes that first element is representative
    elseif X0 isa Number
        d = 1
    else
        throw(ArgumentError("the type of the initial condition, $(typeof(X0)), cannot be handled"))
    end

    if n == d
        return true
    else
        throw_error && throw(ArgumentError("the state-space dimension should match " *
                            "the dimension of the initial state, but they are " *
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
    got_NSTEPS = haskey(kwargs, :NSTEPS)

    if got_tspan && got_T
        throw(ArgumentError("cannot parse the time horizon `T` and the " *
            "time span `tspan` simultaneously; use one or the other"))
    end

    if got_tspan
        tspan =_promote_tspan(kwargs[:tspan])

    elseif got_T
        T = kwargs[:T]
        tspan = (zero(T), T)

    elseif !isempty(args) && (args[1] isa Tuple{Float64, Float64} ||
                              args[1] isa Vector{Float64}         ||
                              args[1] isa Interval || args[1] isa IntervalArithmetic.Interval)
        tspan = _promote_tspan(args[1])

    elseif args[1] isa Float64 # got time horizon as first argument
        tspan = (zero(args[1]), args[1])

    elseif got_NSTEPS
        tspan = (0.0, 0.0) # defined a posteriori

    else
        tspan = (0.0, 0.0)
        # TODO: find better solution such that the error message is used
        #throw(ArgumentError("the time span has not been specified, " *
        #    "but is required for `solve`; you should specify either the time horizon " *
        #    "`T=...` or the time span `tspan=...`"))
    end
    return TimeInterval(tspan[1], tspan[2])
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

const DEFAULT_NSTEPS = 100

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
            opC = GLGM06(δ=δ)
        end
    else
        opC = TMJets()
    end
    return opC
end
