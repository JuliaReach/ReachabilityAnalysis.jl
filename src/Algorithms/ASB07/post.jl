function post(alg::ASB07,
              ivp::IVP{<:AbstractContinuousSystem},
              tspan;
              kwargs...)

    @unpack δ, approx_model, max_order = alg

    # get time horizon from the time span imposing that
    # tspan is of the form (0, T)
    T = _get_T(tspan, check_zero=true, check_positive=true)

    # normalize system to canonical form
    # x' = Ax, x in X
    # x' = Ax + u, x in X, u in U
    ivp_norm = _normalize(ivp)

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)
    Ω0 = initial_state(ivp_discr)
    Φ = state_matrix(ivp_discr)
    X = stateset(ivp_discr)

    Ω0 = _convert_or_overapproximate(Zonotope, Ω0)
    Ω0 = Zonotope(center(Ω0), Matrix(genmat(Ω0)))
    ZT = typeof(Ω0)
    N = eltype(Ω0)

    # preallocate output flowpipe
    NSTEPS = round(Int, T / δ)
    F = Vector{ReachSet{N, ZT}}(undef, NSTEPS)

    if hasinput(ivp_norm)
        error("not implemented")
    end

    reach_homog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X)

    return Flowpipe(F)
end
