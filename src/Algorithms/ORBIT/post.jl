function post(alg::ORBIT{N,VT,AM}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N,VT,AM}
    @unpack δ, approx_model = alg

    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    U = inputset(ivp)

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # preallocate output flowpipe
    ST = Singleton{N,VT}
    F = Vector{ReachSet{N,ST}}(undef, NSTEPS + 1)

    return _post!(U, alg, F, ivp_norm, NSTEPS, Δt0; kwargs...)
end

function _post!(U::AbstractInput, alg::ORBIT, F, ivp, NSTEPS, Δt0; kwargs...)
    U = next_set(U, 1)  # unpack inputs to dispatch on the set type
    return _post!(U, alg, F, ivp, NSTEPS, Δt0; kwargs...)
end

# method for deterministic inputs
function _post!(U::Union{Nothing,AbstractSingleton}, alg::ORBIT, F, ivp,
                NSTEPS, Δt0; kwargs...)
    @unpack δ, approx_model = alg

    # homogenize the initial-value problem
    if get(kwargs, :homogenize, false)
        ivp = homogenize(ivp)
    end

    # discretize system
    ivp_discr = discretize(ivp, δ, approx_model)
    Φ = state_matrix(ivp_discr)
    Ω0 = initial_state(ivp_discr)

    X = stateset(ivp_discr)
    V = inputset(ivp_discr)

    _orbit!(F, Φ, Ω0, V, NSTEPS + 1, δ, X, Δt0)

    return Flowpipe(F)
end

# method for nondeterministic inputs
function _post!(U::LazySet, alg::ORBIT{N,VT}, F, ivp, NSTEPS, Δt0;
                kwargs...) where {N,VT}
    X = stateset(ivp)
    !isa(X, Universe) && error("this algorithm does not " *
                               "support state constraints (also called invariants)")

    @unpack δ, approx_model = alg
    !isa(approx_model, NoBloating) && error("this algorithm does not support " *
                                            "continuous-time analysis")

    A = state_matrix(ivp)
    n = size(A, 1)
    X0 = initial_state(ivp)
    U = next_set(inputset(ivp), 1)

    # discretized matrix
    Φ = _exp(A, δ, approx_model.exp)

    # time interval
    Δt = (zero(N) .. zero(N)) + Δt0

    # start with x0
    x = element(X0)
    @inbounds F[1] = ReachSet(Singleton(x), Δt)
    Δt += δ

    # iterate
    @inbounds for k in 2:(NSTEPS + 1)
        # obtain random input signal and discretize it
        u = sample(U)
        Mu = _Φ₁_u(A, δ, approx_model.exp, approx_model.inv, u, Φ)

        # compute next state
        y = VT(undef, n)
        mul!(y, Φ, x)
        y .+= Mu

        F[k] = ReachSet(Singleton(y), Δt)
        x = y
        Δt += δ
    end

    return Flowpipe(F)
end

# Given a matrix Φ and vectors Ω0 and V, compute the sequence:
#
# Ω0
# Φ*Ω0 + V
# Φ^2*Ω0 + Φ*V + V
# ....
# Φ^(NSTEPS)*Ω0 + Φ^(NSTEPS-1)*V + ... + Φ*V + V
function _orbit!(F, Φ::AbstractMatrix{N}, Ω0, V, NSTEPS, δ, ::Universe, Δt0) where {N}
    # preallocate output sequence
    n = size(Φ, 1)
    VT = vector_type(typeof(Φ))
    out = Vector{VT}(undef, NSTEPS)
    @inbounds for i in 1:NSTEPS
        out[i] = VT(undef, n)
    end

    # compute output sequence
    _orbit!(out, Φ, Ω0, V, NSTEPS)

    # fill reach-set sequence for each time instance 0, δ, 2δ, ...
    Δt = (zero(N) .. zero(N)) + Δt0
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
    @inbounds for i in 2:(NSTEPS - 1)
        mul!(out[i + 1], Φ, out[i])
        out[i + 1] .+= v
    end
    return out
end

# case without input V = 0
function _orbit!(out, Φ::AbstractMatrix{N}, Ω0::Singleton, V::Union{ZeroSet,Nothing},
                 NSTEPS) where {N}
    x = element(Ω0)
    copy!(out[1], x)
    @inbounds for i in 1:(NSTEPS - 1)
        mul!(out[i + 1], Φ, out[i])
    end
    return out
end

# general case
function _orbit!(out, Φ::AbstractMatrix{N}, Ω0::Singleton, V::Singleton, NSTEPS) where {N}
    v = element(V)
    x = element(Ω0)
    copy!(out[1], x)
    @inbounds for i in 1:(NSTEPS - 1)
        mul!(out[i + 1], Φ, out[i])
        out[i + 1] .+= v
    end
    return out
end

function load_krylov_ORBIT()
    return quote

        # case Ω0 = 0 and vector V
        function _orbit_krylov!(A::AbstractMatrix, V::AbstractVector, NSTEPS;
                                hermitian=false, m=min(30, size(A, 1)), tol=1e-7)
            T = eltype(A)
            Ks = KrylovSubspace{T,real(T)}(length(V), m)
            arnoldi!(Ks, A, V; m=m, ishermitian=hermitian, tol=tol)

            out = Vector{typeof(V)}(undef, NSTEPS)
            @inbounds for i in 1:NSTEPS
                out[i] = similar(V)
            end
            out[1] = copy(V)

            @inbounds for i in 1:(NSTEPS - 1)
                expv!(out[i + 1], i * 1.0, Ks)
                out[i + 1] += out[i]
            end
            return out
        end

        # TODO : general case for vectors Ω0 and V
        # ...

    end
end  # quote / load_krylov_ORBIT()
