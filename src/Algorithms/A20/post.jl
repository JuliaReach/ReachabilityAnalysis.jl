# default continuous post

# discrete post
function post(alg::A20{N}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}
    @unpack δ, max_order = alg

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
    static = get(kwargs, :static, false)
    Ω0 = _reconvert(Ω0, static, dim, ngens)
    Φ = _reconvert(Φ, static, dim)

    # preallocate output flowpipe
    #N = eltype(Ω0)
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N,ZT}}(undef, NSTEPS)

    if got_homogeneous

        #=
        # TEMP: static + with preallocation not implemented
        if static == Val(true)
            if alg.preallocate == Val(true)
                @warn "preallocate option is being ignored"
            end
            preallocate = Val(false)
        end
        =#

        reach_homog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, preallocate, Δt0,
                            disjointness_method)
    else
        # TODO: implement preallocate option for this scenario
        U = inputset(ivp)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        U = _convert_or_overapproximate(Zonotope, U)
        reach_inhomog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U, reduction_method, Δt0,
                              disjointness_method)
    end

    return Flowpipe(F)
end
