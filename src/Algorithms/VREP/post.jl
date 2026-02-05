# default continuous post

# discrete post
function post(alg::VREP{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroT, kwargs...) where {N}
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
