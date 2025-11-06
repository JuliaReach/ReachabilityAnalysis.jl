function _get_numsteps(Tsample, δ, ζ, ::DeterministicSwitching)
    αlow = Tsample / δ
    NLOW = ceil(Int, αlow)
    return NLOW, NLOW
end

function _get_numsteps(Tsample, δ, ζ, ::NonDeterministicSwitching)
    ζ⁻ = inf(ζ)
    αlow = (Tsample + ζ⁻) / δ
    NLOW = ceil(Int, αlow)

    ζ⁺ = sup(ζ)
    αhigh = (Tsample + ζ⁺) / δ
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
    αlow = (Tsample + ζ⁻) / δ
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
    if isnothing(alg)
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
        tsp ≠ emptyinterval() ||
            throw(ArgumentError("either the time span `tspan` or the maximum number " *
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
        return ReachSolution(Flowpipe(sol[1:(NLOW - 1)]), alg)
    end

    # preallocate output vector of flowpipes
    N = numtype(alg)
    RT = rsetrep(alg)
    SRT = SubArray{RT,1,Vector{RT},Tuple{UnitRange{Int}},true}
    FT = Flowpipe{N,RT,SRT}
    out = Vector{ShiftedFlowpipe{FT,N}}()
    sizehint!(out, max_jumps + 1)

    push!(out, ShiftedFlowpipe(Flowpipe(view(array(sol), 1:NHIGH)), t0))

    # prepare successor for next jump
    Xend = _transition_successors(sol, NLOW, NHIGH, switching)
    t0 += Tsample + ζ⁻

    # adjust integer bounds for subsequent jumps
    NHIGH += ceil(Int, abs(ζ⁻) / δ)

    @inbounds for k in 2:(max_jumps + 1)

        # apply reset map and instantiate new initial-value problem
        prob = IVP(sys, rmap(Xend))

        # solve next chunk
        sol = post(alg, prob, no_tspan; NSTEPS=NHIGH)

        # store flowpipe until first intersection with the guard
        aux = view(array(sol), 1:NHIGH)

        push!(out, ShiftedFlowpipe(Flowpipe(aux), t0))

        # get successors after discrete jump from Tsample + ζ⁻ .. Tsample + ζ⁺
        Xend = _transition_successors(sol, NLOW, NHIGH, switching)

        # adjust initial time for next jump
        t0 += Tsample
    end

    return ReachSolution(HybridFlowpipe(out), alg)
end

# deterministic switching
@inline function _transition_successors(sol, NLOW, NHIGH, ::DeterministicSwitching)
    # in this scenario, NLOW == NHIGH most of the time
    # sol[NLOW] |> set
    # however, it may happen that NHIGH = NLOW + 1 if Tsample is a multiple of δ
    if NLOW == NHIGH
        # TODO: can we use dispatch and remove this branch?
        set(sol[NLOW])
    else
        set(convexify(Flowpipe(view(array(sol), NLOW:NHIGH))))
    end
end

# non-deterministic switching
@inline function _transition_successors(sol, NLOW, NHIGH, ::NonDeterministicSwitching)
    return set(convexify(Flowpipe(view(array(sol), NLOW:NHIGH))))
end
