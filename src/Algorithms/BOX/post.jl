function post(alg::BOX{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}

    @unpack δ, approx_model, static, dim, recursive = alg

    if haskey(kwargs, :NSTEPS)
        NSTEPS = kwargs[:NSTEPS]
        T = NSTEPS * δ
    else
        # get time horizon from the time span imposing that it is of the form (0, T)
        T = _get_T(tspan, check_zero=true, check_positive=true)
        NSTEPS = ceil(Int, T / δ)
    end

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenize the initial-value problem
    if haskey(kwargs, :homogenize) && kwargs[:homogenize] == true
        ivp_norm = homogenize(ivp_norm)
    end

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)
    Φ = state_matrix(ivp_discr)
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # this algorithm requires Ω0 to be hyperrectangle
    Ω0 = _overapproximate(Ω0, Hyperrectangle)

    # reconvert the set of initial states and state matrix, if needed
    #static = haskey(kwargs, :static) ? kwargs[:static] : alg.stati

    Ω0 = _reconvert(Ω0, static, dim)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    #N = eltype(Ω0)
    HT = typeof(Ω0)
    F = Vector{ReachSet{N, HT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_BOX!(F, Ω0, Φ, NSTEPS, δ, X, recursive, Δt0)
    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expcted input of type `<:LazySet`, but got $(typeof(U))"
        # TODO: can we use support function evaluations for the input set?
        U = overapproximate(U, Hyperrectangle)
        reach_inhomog_BOX!(F, Ω0, Φ, NSTEPS, δ, X, U, recursive, Δt0)
    end

    return Flowpipe(F)
end
