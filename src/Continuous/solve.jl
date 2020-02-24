export solve

# ======================================
# Solve function for continuous systems
# ======================================

function solve(ivp::IVP{<:ACS}, args...; kwargs...)

    S = system(ivp)
    X0 = initial_state(ivp)

    # check dimension
    if statedim(ivp) != dim(X0)
        throw(ArgumentError("the state-space dimension should match the " *
                            "dimension of the initial state, but they are of size " *
                            "$(statedim(problem)) and $(dim(X0)) respectively"))
    end

    tspan = _get_time_span(kwargs...)
    opC = _get_continuous_post(args...)

    _solve_continuous(ivp, opC; tspan=tspan)
end

function _solve_continuous(problem, op; T...)
    # normalize system to canonical form if needed
    # TODO check normalization... probably only applies locally
    problem = IVP(normalize(problem.s), problem.x0)

    # run the continuous-post operator
    sol = post(problem, op; T=T)

    return sol
end

# ==================
# Argument handling
# ==================

@inline function _get_continuous_post(args...)

end

@inline function _get_time_span(; kwargs...)
    got_tspan = haskey(kwargs, :tspan)
    got_T = haskey(kwargs, :T)

    if got_tspan && got_T
        throw(ArgumentError("cannot parse the time horizon `T` and the " *
            "time span `tspan` simultaneously; use one or the other"))
    end

    if got_T
        tspan = kwargs[:tspan]
        tspan = _promote_tspan(tspan)
        # TODO: add functionality for linear systems and remove this restriction
        @assert iszero(tspan[1]) "can only handle initial time t(0) = 0.0"
    elseif got_tspan
        T = kwargs[:T]
        tspan = (zero(T), T)
    else
        throw(ArgumentError("the time span `timespan` has not been specified, " *
            "but is required for `solve`"))
    end
    return tspan
end

@inline _promote_tspan((t1,t2)::Tuple{T,S}) where {T, S} = promote(t1, t2)
@inline _promote_tspan(tspan::Number) = (zero(tspan),tspan)
@inline _promote_tspan(tspan::IA.Interval) = (inf(tspan), sup(tspan))
@inline function _promote_tspan(tspan::AbstractArray)
    if length(tspan) == 2
        return (first(tspan),last(tspan))
    else
        throw(ArgumentError("the length of tspan must be two (and preferably, " *
        "`tspan` should be a tuple, i.e. (0.0,1.0))."))
    end
end

# ===================
# Algorithm defaults
# ===================

"""
    _default_continuous_post(ivp::IVP{ST}) where {ST<:ACS}

Return the default continous post operator for an initial value problem of a
continuous system.

### Input

- `ivp` -- initial value problem

### Output

A continuous post operator with default options.

### Algorithm

- If the system is affine, the algorithm `GLGM06` is used.
- Otherwise, the algorithm `TMJets` is used.
"""
function _default_continuous_post(ivp::IVP{ST}) where {ST<:ACS}
    if isaffine(ivp)
        opC = GLGM06()
    else
        opC = TMJets()
    end
    return opC
end
