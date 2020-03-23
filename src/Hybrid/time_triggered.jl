# ===============================================
# Hybrid systems with time-triggered transitions
# ===============================================

abstract type AbstractHybridAutomatonwithClockedLinearDynamics <: AbstractHybridSystem end
const AHACLD = AbstractHybridAutomatonwithClockedLinearDynamics

"""
    HACLD1{T<:AbstractSystem, ST, MT, N} <: AHACLD

Single-mode hybrid automaton with clocked linear dynamic.

### Fields

- `sys`     -- system
- `X0`      -- initial states
- `rmap`    -- reset map
- `Tsample` -- sampling time
- `ζ`       -- jitter

### Notes

This type is parametric in:

- `T`  -- system type
- `ST` -- type of the initial states
- `MT` -- type of the reset map
- `N`  -- numeric type, for the sampling time and jitter

The following getter functions are available:

- `system`         -- returns the associated system
- `initial_state`  -- returns the set of initial states
- `reset_map`      -- returns the reset map
- `sampling_time`  -- retuns the sampling time
- `jitter`         -- returns the jitter
"""
struct HACLD1{T<:AbstractSystem, ST, MT, N} <: AHACLD
    sys::T
    X0::ST
    rmap::MT
    Tsample::N
    ζ::N
end

MathematicalSystems.system(hs::HACLD1) = hs.sys
MathematicalSystems.initial_state(hs::HACLD1) = hs.x0
reset_map(hs::HACLD1) = hs.rmap
sampling_time(hs::HACLD1) = hs.Tsample
jitter(hs::HACLD1) = hs.ζ

# tstart       Ts-ζ          tend
# [-------------|-------------]
function post(ha::HACLD1,
              max_jumps::Int,
              tspan::TimeInterval,
              alg::AbstractContinuousPost)

    @unpack sys, X0, rmap, Tsample, ζ = ha
    δ = step_size(alg)

    # solve first interval
    prob = IVP(S, X0)
    NLOW = floor(Int, (Tsample - ζ)/δ)
    NHIGH = ceil(Int, (Tsample + ζ)/δ)
    sol = solve(prob, NSTEPS=NHIGH, alg=alg)

    # preallocate output vector of flowpipes
    # TODO: usar un tipo eg. HybridFlowpipe
    FT = typeof(flowpipe(sol))
    out = Vector{FT}()

    @inbounds for k in 1:max_jumps

        # add time interval
        aux = Vector(sol(0 .. Tsample-ζ))

        # porque sol devuelve un view, despues
        # voy a hacer un flowpipe de view
        push!(out, shift(Flowpipe(aux), t0))

        t0 += tstart(aux[end])

        Xend = sol(Tsample-ζ .. Tsample+ζ) |> Vector |> Flowpipe |> Convexify |> set

        sol = solve(IVP(S, reset_map(Xend)), T=Tsample + ζ + δ, alg=alg)

    end

    return out
end
