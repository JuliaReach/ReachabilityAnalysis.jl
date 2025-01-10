# continuous post
function post(alg::INT, ivp::IVP{<:AbstractContinuousSystem}, tspan;
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
function post(alg::INT{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    n = statedim(ivp)
    n == 1 || throw(ArgumentError("this algorithm applies to one-dimensional " *
                                  "systems, but this initial-value problem is $n-dimensional"))

    @unpack δ, approx_model = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)[1, 1]
    Ω0 = initial_state(ivp)
    X = stateset(ivp)

    # this algorithm requires Ω0 to be an interval
    Ω0 = overapproximate(Ω0, Interval)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

    # preallocate output flowpipe
    IT = typeof(Ω0)
    #N = eltype(Ω0)
    F = Vector{ReachSet{N,IT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_INT!(F, Ω0, Φ, NSTEPS, δ, X, Δt0)
    else
        U = inputset(ivp)
        @assert isa(U, LazySet)
        U = overapproximate(U, Interval) # TODO guarantee this on the discretization?
        reach_inhomog_INT!(F, Ω0, Φ, NSTEPS, δ, X, U, Δt0)
    end

    return Flowpipe(F)
end
