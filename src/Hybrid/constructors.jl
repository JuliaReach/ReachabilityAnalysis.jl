# ==========================================
# Convenience hybrid automaton constructors
# ==========================================

import HybridSystems: HybridSystem

# hybrid automaton with one location and a self-loop TODO comparare / use OneStateAutomaton
function HybridSystem(mode::AbstractContinuousSystem, reset_map::AbstractMap)
    automaton = GraphAutomaton(1)
    add_transition!(automaton, 1, 1, 1)
    return HybridSystem(automaton, [mode], [reset_map], [AutonomousSwitching()])
end

# hybrid automaton constructors with default switchings
function HybridSystem(automaton, modes, resetmaps)
    m = nmodes(automaton)
    switchings = fill(AutonomousSwitching(), m)
    return HybridSystem(automaton, modes, resetmaps, switchings)
end

function HybridSystem(; automaton, modes, resetmaps)
    return HybridSystem(automaton, modes, resetmaps)
end

#=
# hybrid automaton with only one mode TODO compare with OneStateAutomaton
struct HA1{ST, IT, TT, GT} <: AbstractHybridSystem
    mode::ST
    invariant::IT
    transition::TT
    guard::GT
end

function HA1(sys::AbstractSystem, tr::AbstractMap)
    # normalize(sys) ? for cases without invariant
    HA1(sys, stateset(sys), tr, stateset(tr))
end
=#

# TODO refactor => MathematicalSystems (?)
HybridSystems.mode(ivp::InitialValueProblem{<:HybridSystem}, i::Integer) = mode(system(ivp), i)

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
                   intervals `[Tsample + inf(ζ), Tsample + sup(ζ)]`; note that
                   the infimum is expected to be negative for most use cases, i.e.
                   when the jitter interval is centered at zero

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

# constructor when jitter is an interval: check whether its diameter is zero or not
# we promote ζ to accept tspan-like inputs (tuples, vectors, etc.)
function HACLD1(sys::T, rmap::MT, Tsample::N, ζ::J) where {T, MT, N, J}
    ζint = _promote_tspan(ζ)
    switching = diam(ζint) > zero(N) ? NonDeterministicSwitching() : DeterministicSwitching()
    return HACLD1(sys, rmap, Tsample, ζint, switching)
end

function _check_dim(ivp::InitialValueProblem{<:HACLD1}; throw_error::Bool=true)
    S = system(ivp)
    X0 = initial_state(ivp)
    _check_dim(S, X0, throw_error=throw_error)
end
