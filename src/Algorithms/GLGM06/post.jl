function post(alg::GLGM06, ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)

    @unpack δ, approx_model, max_order = alg

    # get time horizon from the time span imposing that
    # tspan is of the form (0, T)
    T = _get_T(tspan, check_zero=true, check_positive=true)

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)

    Ω0 = initial_state(ivp_discr)
    Ω0 = _convert_or_overapproximate(Zonotope, Ω0)

    # TEMP keep?
    Ω0 = Zonotope(center(Ω0), Matrix(genmat(Ω0)))
    ZT = Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
    #ZT = typeof(Ω0) ? should be a concretely typed Zonotope

    Φ = state_matrix(ivp_discr)
    N = eltype(Φ)
    X = stateset(ivp_discr)

    # preallocate output flowpipe
    NSTEPS = round(Int, T / δ)
    F = Vector{ReachSet{N, ZT}}(undef, NSTEPS)

    if hasinput(ivp_discr)
        U = inputset(ivp_discr)::LazySet
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U)
    else
        reach_homog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X)
        #reach_homog_inplace!(F, Ω0, Φ, NSTEPS, δ, max_order, X)
    end

    return Flowpipe(F)
end
