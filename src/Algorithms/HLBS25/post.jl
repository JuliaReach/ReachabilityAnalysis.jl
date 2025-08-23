# continuous post
function post(alg::HLBS25, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...)
    δ = alg.δ
    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # discretize system
    ivp_discr = discretize(ivp, δ, alg.approx_model)

    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end

# discrete post
function post(alg::HLBS25{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, max_order, taylor_order, reduction_method, recursive = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)
    Ω0 = initial_state(ivp)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

    # preallocate output flowpipe
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_HLBS25!(F, Ω0, Φ, NSTEPS, δ, max_order, taylor_order, recursive,
                            reduction_method, Δt0)
    else
        error("inhomogeneous algorithm not implemented yet")
    end

    return Flowpipe(F)
end
