function post(alg::ASB07, ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)

    @unpack δ, approx_model, max_order, static, recursive, reduction_method = alg

    if haskey(kwargs, :NSTEPS)
        NSTEPS = kwargs[:NSTEPS]
        T = NSTEPS * δ
    else
        # get time horizon from the time span imposing that it is of the form (0, T)
        T = _get_T(tspan, check_zero=true, check_positive=true)
        NSTEPS = ceil(Int, T / δ)
    end

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
    Ω0 = _reduce_order(Ω0, max_order, reduction_method)

    # reconvert the set of initial states and state matrix, if needed
    #static = haskey(kwargs, :static) ? kwargs[:static] : static # TODO review kwargs vs alg args precedence
    Ω0 = _reconvert(Ω0, static, dim, ngens) # TODO: add dim + generators to the type
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    N = eltype(Ω0)
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N, ZT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X, recursive, reduction_method)
    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expcted input of type `<:LazySet`, but got $(typeof(U))"
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U, reduction_method)
    end

    return Flowpipe(F)
end
