function post(alg::A20{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, max_order = alg

    # TODO move up to main solve function
    if haskey(kwargs, :NSTEPS)
        NSTEPS = kwargs[:NSTEPS]
        T = NSTEPS * δ
    else
        # get time horizon from the time span imposing that it is of the form (0, T)
        T = _get_T(tspan; check_zero=true, check_positive=true)
        NSTEPS = ceil(Int, T / δ)
    end

    # normalize system to canonical form
    # x' = Ax, x in X, or
    # x' = Ax + u, x in X, u in U
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

    # this algorithm requires Ω0 to be a zonotope
    Ω0 = _convert_or_overapproximate(Zonotope, Ω0)
    Ω0 = reduce_order(Ω0, max_order, reduction_method)

    # reconvert the set of initial states and state matrix, if needed
    static = get(kwargs, :static, false)
    Ω0 = _reconvert(Ω0, static, dim, ngens)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    #N = eltype(Ω0)
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)

    if got_homogeneous

        #=
        # TEMP: static + with preallocation not implemented
        if static == Val(true)
            if alg.preallocate == Val(true)
                @warn "preallocate option is being ignored"
            end
            preallocate = Val(false)
        end
        =#

        reach_homog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, preallocate, Δt0,
                            disjointness_method)
    else
        # TODO: implement preallocate option for this scenario
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U, reduction_method, Δt0,
                              disjointness_method)
    end

    return Flowpipe(F)
end
