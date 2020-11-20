# this algorithms assumes that the initial-value problem is one-dimensional
function post(alg::INT{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}

    n = statedim(ivp)
    n == 1 || throw(ArgumentError("this algorithm applies to one-dimensional " *
                    "systems, but this initial-value problem is $n-dimensional"))

    @unpack δ, approx_model = alg

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

    # homogeneize the initial-value problem
    if haskey(kwargs, :homogeneize) && kwargs[:homogeneize] == true
        ivp_norm = homogeneize(ivp_norm)
    end

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)
    Φ = state_matrix(ivp_discr)[1, 1]
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)

    # this algorithm requires Ω0 to be an interval
    Ω0 = overapproximate(Ω0, Interval)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # preallocate output flowpipe
    IT = typeof(Ω0)
    #N = eltype(Ω0)
    F = Vector{ReachSet{N, IT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_INT!(F, Ω0, Φ, NSTEPS, δ, X, Δt0)
    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet)
        U = overapproximate(U, Interval) # TODO guarantee this on the discretization?
        reach_inhomog_INT!(F, Ω0, Φ, NSTEPS, δ, X, U, Δt0)
    end

    return Flowpipe(F)
end
