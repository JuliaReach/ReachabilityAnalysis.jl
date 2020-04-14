# ===============================================
# Hybrid systems with time-triggered transitions
# ===============================================

abstract type AbstractHybridAutomatonwithClockedLinearDynamics <: AbstractHybridSystem end
const AHACLD = AbstractHybridAutomatonwithClockedLinearDynamics

import MathematicalSystems: system, statedim, initial_state

"""
    HACLD1{T<:AbstractSystem, MT, N} <: AHACLD

Single-mode hybrid automaton with clocked linear dynamic.

### Fields

- `sys`     -- system
- `rmap`    -- reset map
- `Tsample` -- sampling time
- `ζ`       -- jitter

### Notes

This type is parametric in:

- `T`  -- system type
- `MT` -- type of the reset map
- `N`  -- numeric type, applies to the sampling time and jitter

The following getter functions are available:

- `initial_state`  -- initial state of the continuous mode
- `jitter`         -- returns the jitter
- `reset_map`      -- returns the reset map
- `sampling_time`  -- retuns the sampling time
- `statedim`       -- dimension of the state-space
- `system`         -- returns the continuous mode
"""
struct HACLD1{T<:AbstractSystem, MT, N} <: AHACLD
    sys::T
    rmap::MT
    Tsample::N
    ζ::N
end

# default constructor without jitter
function HACLD1(sys::T, rmap::MT, Tsample::N) where {T, MT, N}
    return HACLD1(sys, rmap, Tsample, zero(N))
end

initial_state(hs::HACLD1) = initial_state(hs.sys)
jitter(hs::HACLD1) = hs.ζ
reset_map(hs::HACLD1) = hs.rmap
sampling_time(hs::HACLD1) = hs.Tsample
statedim(hs::HACLD1) = statedim(hs.sys)
system(hs::HACLD1) = hs.sys

# tstart       Ts-ζ          tend
# [-------------|-------------]
function post(alg::AbstractContinuousPost, ivp::IVP{<:HACLD1}, tspan; kwargs...)

    X0 = initial_state(ivp)
    ha = system(ivp)
    @unpack sys, rmap, Tsample, ζ = ha
    t0 = tstart(tspan)
    δ = step_size(alg)

    if haskey(kwargs, :max_jumps)
        max_jumps = kwargs[:max_jumps]
    else
        # define max_jumps using the time horizon tspan
        T = tend(tspan)
        # TODO: double check
        α = T / (Tsample - ζ)
        if α <= 0
            error("inconsistent choice of parameters: T / (Tsample - ζ) = $α, " *
                  "but it should be positive")
        end
        max_jumps = ceil(Int, α)
    end

    # solve first interval
    prob = IVP(sys, X0)
    αlow = (Tsample - ζ)/δ
    NLOW = ceil(Int, αlow)
    if NLOW == 0
        error("inconsistent choice of parameters: (Tsample - ζ)/δ = $αlow " *
              "but it should be positive")
    end
    αhigh = (Tsample + ζ)/δ
    NHIGH = ceil(Int, αhigh)
    println(NHIGH)
    sol = solve(prob, NSTEPS=NHIGH, alg=alg; kwargs...)

    # preallocate output vector of flowpipes
    FT = typeof(flowpipe(sol))
    out = Vector{FT}()
    sizehint!(out, max_jumps+1)

    @inbounds for k in 1:max_jumps+1

        # add time interval
        aux = sol[1:NLOW] # Vector(sol(0 .. Tsample-ζ))

        push!(out, shift(Flowpipe(aux), t0))

        t0 += tstart(aux[end])

        # Xend = sol(Tsample-ζ .. Tsample+ζ) |> Vector |> Flowpipe |> Convexify |> set
        Xend = sol[NLOW:NHIGH] |> Vector |> Flowpipe |> Convexify |> set

        prob = IVP(sys, rmap(Xend))
        sol = solve(prob, NSTEPS=NHIGH, alg=alg; kwargs...)
    end

    return HybridFlowpipe(out)
end

#=
function reach_HACLD1_no_jitter()

end

function reach_HACLD1_with_jitter()

end
=#
