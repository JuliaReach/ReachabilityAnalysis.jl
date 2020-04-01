# continuous post for GLGM06 using Zonotope set representation
function post(alg::GLGM06, ivp::IVP{<:AbstractContinuousSystem}, tspan; kwargs...)

    @unpack δ, approx_model, max_order = alg

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
    Φ = state_matrix(ivp_discr)
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # this algorithm requires Ω0 to be a zonotope
    Ω0 = _convert_or_overapproximate(Zonotope, Ω0)
    Ω0 = reduce_order(Ω0, max_order)

    # reconvert the set of initial states, if needed
    force_static = haskey(kwargs, :force_static) ? kwargs[:force_static] : false
    Ω0 = _reconvert(Ω0, Val(force_static))

    if force_static
        n = size(Φ, 1)
        N = eltype(Ω0)
        Φ = SMatrix{n, n, N, n*n}(Φ)
    end

    # preallocate output flowpipe
    N = eltype(Ω0)
    ZT = typeof(Ω0)
    F = Vector{ReachSet{N, ZT}}(undef, NSTEPS)

    if got_homogeneous
        reach_homog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X)
    else
        U = inputset(ivp_discr)
        if isa(U, LazySet)
            U = _convert_or_overapproximate(Zonotope, U)
            reach_inhomog_GLGM06!(F, Ω0, Φ, NSTEPS, δ, max_order, X, U)
        else
            error("inputs of type $(typeof(U)) cannot be handled yet")
        end
    end

    return Flowpipe(F)
end

function _reconvert(Ω0::Zonotope{N, VN, <:Diagonal}, static::Val{false}) where {N, VN}
    c0 = Vector(Ω0.center)
    G0 = Matrix(Ω0.generators)
    Ω0 = Zonotope(c0, G0)
end

function _reconvert(Ω0::Zonotope{N, VN, <:Diagonal}, static::Val{true}) where {N, VN}
    c0 = Ω0.center
    G0 = Ω0.generators
    n = size(G0, 1) # dimension
    p = size(G0, 2) # number of generators
#   Φ = SMatrix{n, n, N, n*n}(Φ)
    c0_st = SVector{n, N}(c0)
    G0_st = SMatrix{n, p, N, n*p}(G0)
    Ω0 = Zonotope(c0_st, G0_st)
end
