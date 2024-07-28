function post(alg::BOX{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, static, dim, recursive = alg

    NSTEPS = _get_nsteps(kwargs, δ, tspan)

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenize the initial-value problem
    if get(kwargs, :homogenize, false)
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
    Ω0 = _reconvert(Ω0, static, dim)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    #N = eltype(Ω0)
    HT = typeof(Ω0)
    F = Vector{ReachSet{N,HT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_BOX!(F, Ω0, Φ, NSTEPS, δ, X, recursive, Δt0)
    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        # TODO: can we use support function evaluations for the input set?
        U = overapproximate(U, Hyperrectangle)
        reach_inhomog_BOX!(F, Ω0, Φ, NSTEPS, δ, X, U, recursive, Δt0)
    end

    return Flowpipe(F)
end
