# continuous post
function post(alg::ASB07, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...)
    δ = alg.δ

    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenization for interval matrices is not available yet

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, alg.approx_model)

    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end

# discrete post
function post(alg::ASB07{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, max_order, reduction_method, static, recursive, dim, ngens = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)
    Ω0 = initial_state(ivp)
    X = stateset(ivp)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

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
        U = inputset(ivp)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_ASB07!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U, recursive, reduction_method, Δt0)
    end

    return Flowpipe(F)
end
