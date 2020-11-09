using LazySets.Arrays: _vector_type

function post(alg::ORBIT{N, VT, AM}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N, VT, AM}

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
    Φ = state_matrix(ivp_discr)
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)
    V = inputset(ivp_discr)

    # preallocate output flowpipe
    ST = Singleton{N, VT}
    F = Vector{ReachSet{N, ST}}(undef, NSTEPS)

    _orbit!(F, Φ, Ω0, V, NSTEPS, δ, X, Δt0)

    return Flowpipe(F)
end

# Given a matrix Φ and vectors Ω0 and V, compute the sequence:
#
# Ω0
# Φ*Ω0 + V
# Φ^2*Ω0 + Φ*V + V
# ....
# Φ^(NSTEPS-1)*Ω0  + Φ^(NSTEPS-2)*V + ... + Φ*V + V
function _orbit!(F, Φ::AbstractMatrix{N}, Ω0, V, NSTEPS, δ, X::Universe, Δt0) where {N}
    # preallocate output sequence
    n = size(Φ, 1)
    VT = _vector_type(typeof(Φ))
    out = Vector{VT}(undef, NSTEPS)
    @inbounds for i in 1:NSTEPS
        out[i] = VT(undef, n)
    end

    # compute output sequence
    _orbit!(out, Φ, Ω0, V, NSTEPS)

    # fill singleton reach-set sequence
    Δt = (zero(N) .. δ) + Δt0
    for k in 1:NSTEPS
        xk = Singleton(out[k])
        F[k] = ReachSet(xk, Δt)
        Δt += δ
    end
    return out
end

# case with zero homogeneous solution Ω0 = 0
function _orbit!(out, Φ::AbstractMatrix{N}, Ω0::ZeroSet, V::Singleton, NSTEPS) where {N}
    v = element(V)
    out[1] = zeros(N, n)
    copy!(out[2], v)
    @inbounds for i in 2:NSTEPS-1
        mul!(out[i+1], Φ, out[i])
        out[i+1] .+= v
    end
    return out
end

# case without input V = 0
function _orbit!(out, Φ::AbstractMatrix{N}, Ω0::Singleton, V::Union{ZeroSet, Nothing}, NSTEPS) where {N}
    x = element(Ω0)
    copy!(out[1], x)
    @inbounds for i in 1:NSTEPS-1
        mul!(out[i+1], Φ, out[i])
    end
    return out
end

# general case
function _orbit!(out, Φ::AbstractMatrix{N}, Ω0::Singleton, V::Singleton, NSTEPS) where {N}
    v = element(V)
    x = element(Ω0)
    copy!(out[1], x)
    @inbounds for i in 1:NSTEPS-1
        mul!(out[i+1], Φ, out[i])
        out[i+1] .+= v
    end
    return out
end
