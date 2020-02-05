export solve

# ===================
# Algorithm defaults
# ===================

"""
    default_continuous_post(ivp::IVP{ST}) where {ST<:AbstractContinuousSystem}

Return the default continous post operator for an initial value problem of a
continuous system.

### Input

- `ivp` -- initial value problem

### Output

A continuous post operator with default options.

### Algorithm

If the system is affine, the algorithm `GLGM06` is returned. Otherwise, the algorithm
`TMJets` is used.
"""
function default_continuous_post(ivp::IVP{ST}) where {ST<:AbstractContinuousSystem}
    if isaffine(ivp)
        opC = GLGM06()
    else
        opC = TMJets()
    end
    return opC
end

# main solve function for continuous systems
function solve(ivp::InitialValueProblem{<:AbstractContinuousSystem}, args...; kwargs...)

    _check_dimension(ivp)
    tspan = _get_time_span(kwargs...)
    opC = _get_continuous_post(args...)
    _solve_continuous(ivp, opC; tspan=tspan)
end

@inline function _check_dimension(ivp::InitialValueProblem)
    X0 = initial_state(ivp)
    if statedim(ivp) == dim(X0)
        return true
    else
        throw(ArgumentError("the state-space dimension should match the " *
                            "dimension of the initial state, but they are of size " *
                            "$(statedim(problem)) and $(dim(X0)) respectively"))
    end
end

# either T or tspan should be specified
@inline function _get_time_span(; kwargs...)
    tspan_specified = haskey(kwargs, :tspan)
    time_horizon__specified = haskey(kwargs, :T)
    if tspan_specified && time_horizon_specified
        throw(ArgumentError("cannot specify both the time horizon `T` and the time span `tspan` at the same time"))
    elseif tspan_specified
        tspan = kwargs[:tspan]
        @assert iszero(tspan[1])
    elseif time_horizon_specified
        T = kwargs[:T]
        tspan = (zero(T), T)
    else
        throw(ArgumentError("the time span `timespan` has not been specified, but is required for `solve`"))
    end
end

@inline function _get_continuous_post_operator(args...)

end

function _solve_continuous(problem, op; T...)
    # normalize system to canonical form if needed
    # TODO check normalization... probably only applies locally
    problem = IVP(normalize(problem.s), problem.x0)

    # run the continuous-post operator
    sol = post(problem, op; T=T)

    return sol
end

#init(::ProblemType, args...; kwargs...) :: SolverType
#solve!(::SolverType) :: SolutionType
