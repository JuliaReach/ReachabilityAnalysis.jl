# this algorithms assumes that the initial-value problem
# is one-dimensional
function post(alg::INT20, ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)

    @assert statedim(ivp) == 1 "this algrithm only applies to one-dimensional systems, " *
                               "but this system is $(statedim(ivp))-dimensional"

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

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)
    Φ = state_matrix(ivp_discr)[1, 1]  # can use first(...) in Julia v1.4
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)

    # this algorithm requires Ω0 to be an interval
    Ω0 = overapproximate(Ω0, Interval)

    # true <=> there is no input, i.e. the system is of the form
    # x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # preallocate output flowpipe
    IT = typeof(Ω0)
    N = eltype(Ω0)
    F = Vector{ReachSet{N, IT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_INT20!(F, Ω0, Φ, NSTEPS, δ, X)
    else
        error("not implemented yet")
        U = inputset(ivp_discr)
        if isa(U, LazySet)
            U = overapproximate(U, Interval)
            reach_inhomog_INT20!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U)
        else
            error("inputs of type $(typeof(U)) cannot be handled yet")
        end
    end

    return Flowpipe(F)
end
