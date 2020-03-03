# TODO: dispatch on linear system only.. ?
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
    ZT = Zonotope{Float64, Vector{Float64}, Matrix{Float64}} # TODO: typeof(Ω0) # should be a concretely typed Zonotope
    F = Vector{ReachSet{N, ZT}}(undef, NSTEPS)
    if hasinput(ivp)
        U = inputset(ivp_discr)::LazySet
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog!(F, Ω0, Φ, NSTEPS, δ, alg.max_order, X, U)
        # error("time-varying input sets not implemented yet")
    else
        reach_homog!(F, Ω0, Φ, NSTEPS, δ, alg.max_order, X)
    end

    return Flowpipe(F)
end
