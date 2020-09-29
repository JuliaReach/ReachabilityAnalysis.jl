# this algorithms assumes that the initial-value problem is one-dimensional
# FIXME requires special scalar quadratic system type
function reach_homog_QINT(alg::QINT{N}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}

    n = statedim(ivp)
    n == 1 || throw(ArgumentError("this algorithm applies to one-dimensional " *
                    "systems, but this initial-value problem is $n-dimensional"))

    @unpack Δ, δ, θ, maxiter, approx_model = alg

    # normalize system to canonical form FIXME
    ivp_norm = _normalize(ivp)
    Ω0 = initial_state(ivp_norm)

    # discretize system
    #ivp_discr = discretize(ivp_norm, δ, approx_model)
    #Φ = state_matrix(ivp_discr)[1, 1]
    #Ω0 = initial_state(ivp_discr)
    #X = stateset(ivp_discr)

    # this algorithm requires Ω0 to be an interval
    #Ω0 = overapproximate(Ω0, Interval)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # preallocate output flowpipe
    #IT = typeof(Ω0)
    #N = eltype(Ω0)
    #F = Vector{ReachSet{N, IT}}(undef, NSTEPS)

    if got_homogeneous
        # FIXME consider Δt0
        F = reach_homog_QINT(a=a, b=b, c=c, X0=Ω0, T=T, Δ=Δ, δ=δ, θ=θ, maxiter=maxiter)
    else
        error("not implemented yet")
        #U = inputset(ivp_discr)
        #@assert isa(U, LazySet)
        #U = overapproximate(U, Interval) # TODO guarantee this on the discretization?
        #reach_inhomog_INT!(F, Ω0, Φ, NSTEPS, δ, X, U, Δt0)
    end

    return F
end
