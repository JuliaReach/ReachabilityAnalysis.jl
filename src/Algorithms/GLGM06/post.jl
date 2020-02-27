# return the time horizon given a time span
function _get_T(tspan::Tuple{Float64, Float64})
    t0 = tspan[1]
    # TODO: add this functionality, see #21
    @assert iszero(t0) "this algorithm can only handle zero initial time"

    T = tspan[2]
    @assert T > 0 "the time horizon should be positive"

    return T
end

# The canonical form is:
#     - If the system doesn't have input, a constrained linear continous system (CLCS)
#       x' = Ax, x ∈ X
#     - If the system has an input, a CLCCS, x' = Ax + u, x ∈ X, u ∈ U
# If the original system is unconstrained, the constraint set X is the universal set.
function _normalize(ivp)

    # initial states normalization
    X0 = initial_state(ivp)
    if X0 isa AbstractVector
        X0_norm = Singleton(X0)
    elseif X0 isa IA.Interval
        X0_norm = convert(Interval, X0)
    elseif X0 isa IA.IntervalBox
        X0_norm = convert(Hyperrectangle, X0)
    else
        X0_norm = X0
    end

    # system's normalization
    S = system(ivp)
    S_norm = normalize(S)

    if S_norm === S && X0_norm === X0
        ivp_norm = ivp
    else
        ivp_norm = IVP(S_norm, X0_norm)
    end

    return ivp_norm
end

# TODO: refactor
hasinput(S::AbstractSystem) = inputdim(S) > 0
isconstantinput(::ConstantInput) = true
isconstantinput(::VaryingInput) = false
isconstantinput(::LazySet) = true

function post(alg::GLGM06, ivp::IVP{<:AbstractContinuousSystem}, tspan, args...; kwargs...)

    # get time horizon
    T = _get_T(tspan)

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # discretize system
    δ = step_size(alg)
    ivp_discr = discretize(ivp_norm, δ, alg.approximation_model)
    Ω0 = initial_state(ivp_discr)
    Ω0 = _convert_or_overapproximate(Zonotope, Ω0)
    Φ = state_matrix(ivp_discr)
    N = eltype(Φ)

    # flowpipe computation
    NSTEPS = round(Int, T / δ)
    F = Vector{ReachSet{N, Zonotope{N}}}(undef, NSTEPS)
    if hasinput(ivp)
        U = inputset(ivp_discr)::LazySet
        reach_inhomog!(F, Ω0, Φ, NSTEPS, δ, alg.max_order, U)
        # error("time-varying input sets not implemented yet")
    else
        reach_homog!(F, Ω0, Φ, NSTEPS, δ, alg.max_order)
    end

    return Flowpipe(F)
end
