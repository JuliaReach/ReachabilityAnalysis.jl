export solve

# ======================================
# Solve function for continuous systems
# ======================================

function solve(ivp::IVP{<:AbstractContinuousSystem}, args...; kwargs...)
    _check_dim(ivp)
    tspan = _get_tspan(; kwargs...)
    cpost = _get_cpost(ivp, tspan, args...; kwargs...)

    # run the continuous-post operator
    post(cpost, ivp, tspan, args...; kwargs...)
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

function _check_dim(ivp)
    S = system(ivp)
    X0 = initial_state(ivp)
    if statedim(S) != dim(X0)
        throw(ArgumentError("the state-space dimension should match the " *
                            "dimension of the initial state, but they are of size " *
                            "$(statedim(S)) and $(dim(X0)) respectively"))
    else
        return true
    end
end

@inline _promote_tspan((t1, t2)::Tuple{T, S}) where {T, S} = promote(t1, t2)
@inline _promote_tspan(tspan::Number) = (zero(tspan),tspan)
@inline _promote_tspan(tspan::IA.Interval) = (inf(tspan), sup(tspan))
@inline function _promote_tspan(tspan::AbstractArray)
    if length(tspan) == 2
        return (first(tspan), last(tspan))
    else
        throw(ArgumentError("the length of tspan must be two (and preferably, " *
                            "`tspan` should be a tuple, i.e. (0.0,1.0)), but " *
                            "it is of length $(length(tspan))"))
    end
end

function _get_tspan(; kwargs...)
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
    return tspan
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
            δ = (tspan[2] - tspan[1]) / N
        else
            δ = (tspan[2] - tspan[1]) / DEFAULT_NSTEPS
        end
        opC = GLGM06(δ=δ)
    else
        opC = TMJets()
    end
    return opC
end
