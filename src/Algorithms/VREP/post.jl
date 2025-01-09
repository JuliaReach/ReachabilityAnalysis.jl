# continuous post
function post(alg::VREP, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...)
    δ = alg.δ

    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenize the initial-value problem
    if get(kwargs, :homogenize, false)
        ivp_norm = homogenize(ivp_norm)
    end

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, alg.approx_model)

    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end

# discrete post
function post(alg::VREP{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, static, dim = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)

    # the initial state should be expressed in v-rep
    Ω0 = initial_state(ivp)
    if !(Ω0 isa VPolygon || Ω0 isa VPolytope)
        n = size(Φ, 1)
        Ω0 = n == 2 ? convert(VPolygon, Ω0) : convert(VPolytope, Ω0)
    end

    X = stateset(ivp)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

    # reconvert the set of initial states and state matrix, if needed
    Ω0 = _reconvert(Ω0, static, dim)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    VT = typeof(Ω0)
    F = Vector{ReachSet{N,VT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_VREP!(F, Ω0, Φ, NSTEPS, δ, X, Δt0)
    else
        error("not implemented")
    end

    return Flowpipe(F)
end
