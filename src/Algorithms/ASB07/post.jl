function post(alg::ASB07{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, max_order, reduction_method, static, recursive, dim, ngens = alg

    NSTEPS = _get_nsteps(kwargs, δ, tspan)

    # normalize system to canonical form
    # x' = Ax, x in X, or
    # x' = Ax + u, x in X, u in U
    ivp_norm = _normalize(ivp)

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
    Ω0 = _reconvert(Ω0, static, dim, ngens)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    #N = eltype(Ω0)
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X, recursive, reduction_method, Δt0)
    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U, recursive, reduction_method, Δt0)
    end

    return Flowpipe(F)
end
