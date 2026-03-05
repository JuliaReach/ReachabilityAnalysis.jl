# continuous post
function post(alg::HLBS25, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroT, kwargs...)
    δ = alg.δ
    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # discretize system
    ivp_discr = discretize(ivp, δ, alg.approx_model)

    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end

# discrete post
function post(alg::HLBS25{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, approx_model, taylor_order, max_order, max_order_zono, reduction_method, recursive, idg = alg

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

    if got_homogeneous
        # homogeneous branch keeps a fixed sparse-polynomial-zonotope representation
        ZT = typeof(Ω0)
        F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)
        reach_homog_HLBS25!(F, Ω0, Φ, NSTEPS, δ, taylor_order, recursive, max_order, max_order_zono,
                            reduction_method, Δt0, idg)
    else
        # concretize the initial exact-sum container to obtain the concrete set type
        Ω0c = concretize(Ω0)
        ZT = typeof(Ω0c)
        F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)
        B = input_matrix(ivp)
        Φ_norm = norm(Φ, Inf)
        U = inputset(ivp)
        reach_inhomog_HLBS25!(F, Ω0, Φ, B, U, NSTEPS, δ, taylor_order, Φ_norm, recursive, max_order,
                              max_order_zono, reduction_method, Δt0, idg)
    end

    return Flowpipe(F)
end
