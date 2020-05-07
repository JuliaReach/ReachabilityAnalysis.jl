# ===============================================
# Hybrid systems with time-triggered transitions
# ===============================================

import MathematicalSystems: system, statedim, initial_state

# particular HA where transitions are only time-triggered
abstract type AbstractHybridAutomatonwithClockedLinearDynamics <: AbstractHybridSystem end

# type alias
const AHACLD = AbstractHybridAutomatonwithClockedLinearDynamics

# abstarct supertype for the types of switching between discrete modes that are
# handled by hybrid systems with time-triggeerd transitions
abstract type AbstractSwitching end

# deterministic switchings occur at prescribed time stamps, i.e. there is no jitter
struct DeterministicSwitching <: AbstractSwitching end

#const deterministic_switching = DeterministicSwitching()

# non-deterministic switchings may occur at any time between a prescribed interval,
# i.e. there jitter is represented by an interval
struct NonDeterministicSwitching <: AbstractSwitching end

#const nondeterministic_switching() = DeterministicSwitching()

"""
    HACLD1{T<:AbstractSystem, MT, N, J} <: AHACLD

Single-mode hybrid automaton with clocked linear dynamics.

### Fields

- `sys`       -- system
- `rmap`      -- reset map
- `Tsample`   -- sampling time
- `ζ`         -- jitter
- `switching` -- (optional) value that

### Notes

This type is parametric in:

- `T`  -- system type
- `MT` -- type of the reset map
- `N`  -- numeric type, applies to the sampling time and jitter
- `J`  -- type associated to the jitter

The type associated to the jitter, `J`, can be one of the following:

- `Missing`     -- no jitter, i.e. switchings are deterministic
- `Number`      -- symetric jitter, i.e. non-deterministic switchings occur in the
                   intervals `[Tsample - ζ, Tsample + ζ]`
- `IA.Interval` -- nonsymetric jitter, i.e. non-deterministic switchings occur in the
                   intervals `[Tsample - inf(ζ), Tsample + sup(ζ)]`

The following getter functions are available:

- `initial_state`  -- initial state of the continuous mode
- `jitter`         -- return the jitter
- `reset_map`      -- return the reset map
- `sampling_time`  -- return the sampling time
- `statedim`       -- dimension of the state-space
- `system`         -- return the continuous mode
- `switching`      -- return the type of switching

Non-deterministic switching:

 tstart       Ts-ζ⁻         tend
 [-------------|-------------]

In the following, suppose that the continuous post-operator has fixed step-size
`δ > 0`. If `F` denotes the flowpipe, then

  F[1]    F[2]    F[3]    F[4]    F[5]        F[k]
[------][------][------][------][------]  ⋅ ⋅ ⋅ ⋅ ⋅ ⋅ ⋅ ⋅  [------]

 `R = array(F)` denote the array of reach-sets time-span for reach reach-set is of the form:

Similarly we compute tstart and tend for the supremum part Ts-ζ⁺
"""
struct HACLD1{T<:AbstractSystem, MT, N, J, S<:AbstractSwitching} <: AHACLD
    sys::T
    rmap::MT
    Tsample::N
    ζ::J
    switching::S
end

# getter functions
system(hs::HACLD1) = hs.sys
reset_map(hs::HACLD1) = hs.rmap
sampling_time(hs::HACLD1) = hs.Tsample
statedim(hs::HACLD1) = statedim(hs.sys)
initial_state(hs::HACLD1) = initial_state(hs.sys)
jitter(hs::HACLD1{T, MT, N, Missing}) where {T, MT, N} = TimeInterval(zero(N), zero(N))
jitter(hs::HACLD1{T, MT, N, J}) where {T, MT, N, J}  = hs.ζ
# TODO cleanup
#jitter(hs::HACLD1{T, MT, N, J}) where {T, MT, N, J<:Real} = TimeInterval(-ζ, ζ)
#jitter(hs::HACLD1{T, MT, N, J}) where {T, MT, N, J}  = _promote_tspan(hs.ζ)
switching(::HACLD1{T, MT, N, J, S}) where {T, MT, N, J, S} = S

# default constructor without jitter
function HACLD1(sys::T, rmap::MT, Tsample::N) where {T, MT, N}
    return HACLD1(sys, rmap, Tsample, missing, DeterministicSwitching())
end

# constructor when jitter is a real number
function HACLD1(sys::T, rmap::MT, Tsample::N, ζ::J) where {T, MT, N, J<:Real}
    if iszero(ζ)
        switching = DeterministicSwitching()
        ζint = TimeInterval(ζ, ζ)
    else
        switching = NonDeterministicSwitching()
        ζint = TimeInterval(-ζ, ζ)
    end
    return HACLD1(sys, rmap, Tsample, ζint, switching)
end

# should be treated separately because IA.Interval <: Number
function HACLD1(sys::T, rmap::MT, Tsample::N, ζ::J) where {T, MT, N, J<:IA.Interval}
    return HACLD1(sys, rmap, Tsample, ζ, NonDeterministicSwitching())
end

# constructor when jitter is an inteval: check whether its diameter is zero or not
# we promote ζ to accept tspan-like inputs (tuples, vectors, etc.)
function HACLD1(sys::T, rmap::MT, Tsample::N, ζ::J) where {T, MT, N, J}
    ζint = _promote_tspan(ζ)
    switching = diam(ζint) > zero(N) ? NonDeterministicSwitching() : DeterministicSwitching()
    return HACLD1(sys, rmap, Tsample, ζint, switching)
end

function _get_numsteps(Tsample, δ, ζ, ::DeterministicSwitching)
    αlow = Tsample / δ
    NLOW = ceil(Int, αlow)
    return NLOW, NLOW
end

function _get_numsteps(Tsample, δ, ζ, ::NonDeterministicSwitching)
    ζ⁻ = inf(ζ)
    αlow = (Tsample + ζ⁻)/δ
    NLOW = ceil(Int, αlow)

    ζ⁺ = sup(ζ)
    αhigh = (Tsample + ζ⁺)/δ
    NHIGH = ceil(Int, αhigh)

    return NLOW, NHIGH
end

function _get_max_jumps(tspan, Tsample, δ, ζ, ::DeterministicSwitching)
    NLOW = ceil(Tsample / δ)
    T = tend(tspan)
    aux = (T - (NLOW - 1) * δ) / (NLOW * δ)
    max_jumps = ceil(Int, aux)
    @assert max_jumps >= 0 "inconsistent choice of parameters, `max_jumps` should " *
                           "be positive but it is $max_jumps"

    return max_jumps
end

function _get_max_jumps(tspan, Tsample, δ, ζ, ::NonDeterministicSwitching)
    ζ⁻ = inf(ζ)
    αlow = (Tsample + ζ⁻)/δ
    NLOW = ceil(Int, αlow)
    T = tend(tspan)
    aux = (T - (NLOW - 1) * δ) / (NLOW * δ)
    max_jumps = ceil(Int, aux)
    @assert max_jumps >= 0 "inconsistent choice of parameters, `max_jumps` should " *
                           "be positive but it is $max_jumps"

    return max_jumps
end

const no_tspan = emptyinterval()

function solve(ivp::IVP{<:HACLD1}, args...; kwargs...)

    # preliminary checks
    _check_dim(ivp)

    # get time span (or the emptyset if NSTEPS was specified)
    tsp = _get_tspan(args...; kwargs...)

    # get the continuous post or find a default one
    alg = _get_cpost(ivp, args...; kwargs...)
    if alg == nothing
        alg = _default_cpost(ivp, tsp; kwargs...)
    end

    X0 = initial_state(ivp)
    ha = system(ivp)
    @unpack sys, rmap, Tsample, ζ, switching = ha
    t0 = isempty(tsp) ? 0.0 : tstart(tsp)
    δ = step_size(alg)

    # get maximum number of jumps
    if haskey(kwargs, :max_jumps)
        max_jumps = kwargs[:max_jumps]
    else
        tsp ≠ emptyinterval() || throw(ArgumentError("either the time span `tspan` or the maximum number " *
                                                       "of jumps `max_jumps` should be specified"))
        max_jumps = _get_max_jumps(tsp, Tsample, δ, ζ, switching)
    end

    # number of steps for the continuous post
    prob = IVP(sys, X0)
    NLOW, NHIGH = _get_numsteps(Tsample, δ, ζ, switching)
    ζint = jitter(ha)
    ζ⁻ = inf(ζint)
    ζ⁺ = sup(ζint)
    @assert NLOW > 0 throw(ArgumentError("inconsistent choice of parameters"))

    # solve first flowpipe
    sol = post(alg, prob, no_tspan; NSTEPS=NHIGH)

    if max_jumps == 0
        return ReachSolution(Flowpipe(sol[1:NLOW-1]), alg)
    end

    # preallocate output vector of flowpipes
    N = numtype(alg)
    RT = rsetrep(alg)
    SRT = SubArray{RT, 1, Vector{RT}, Tuple{UnitRange{Int}}, true}
    FT = Flowpipe{N, RT, SRT}
    out = Vector{ShiftedFlowpipe{FT, N}}()
    sizehint!(out, max_jumps+1)

    push!(out, ShiftedFlowpipe(Flowpipe(view(array(sol), 1:NHIGH)), t0))

    # prepare successor for next jump
    Φ = exp(state_matrix(sys) * Tsample)
    Xend = _transition_successors(ivp, X0, Φ)
    t0 += Tsample + ζ⁻

    # adjust integer bounds for subsequent jumps
    NHIGH += ceil(Int, abs(ζ⁻) / δ)

    @inbounds for k in 2:max_jumps+1

        # apply reset map and instantiate new initial-value problem
        X0_new = rmap(Xend)
        prob = IVP(sys, X0_new)

        # solve next chunk
        sol = post(alg, prob, no_tspan; NSTEPS=NHIGH)

        # store flowpipe until first intersection with the guard
        aux = view(array(sol), 1:NHIGH)

        push!(out, ShiftedFlowpipe(Flowpipe(aux), t0))

        # get successors after discrete jump from Tsample + ζ⁻ .. Tsample + ζ⁺
        Xend = _transition_successors(ivp, X0_new, Φ)

        # adjust initial time for next jump
        t0 += Tsample
    end

    return ReachSolution(HybridFlowpipe(out), alg)
end

# deterministic switching
@inline function _transition_successors(ivp, X0, Φ)
    return linear_map(Φ, X0)
end

# non-deterministic switching
@inline function _transition_successors(sol, NLOW, NHIGH, ::NonDeterministicSwitching)
    view(array(sol), NLOW:NHIGH) |> Flowpipe |> Convexify |> set
end
