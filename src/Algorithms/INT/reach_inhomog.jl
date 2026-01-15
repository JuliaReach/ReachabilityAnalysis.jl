# ===================
# Inhomogeneous case
# ===================

# inhomogeneous case; no invariant; recursive implementation
function reach_inhomog_INT!(F::Vector{ReachSet{N,Interval{N}}},
                            Ω0::Interval{N},
                            Φ::N,
                            NSTEPS::Integer,
                            δ::Float64,
                            ::Universe,
                            U::Interval,
                            Δt0::TimeInterval) where {N}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Iₖ = Interval(Φ * set(F[k - 1]).dat + U.dat)
        Δt += δ
        F[k] = ReachSet(Iₖ, Δt)
        k += 1
    end
    return F
end

# inhomogeneous case; with invariant; recursive implementation
function reach_inhomog_INT!(F::Vector{ReachSet{N,Interval{N}}},
                            Ω0::Interval{N},
                            Φ::N,
                            NSTEPS::Integer,
                            δ::Float64,
                            X::LazySet,
                            U::Interval,
                            Δt0::TimeInterval) where {N}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Iₖ = Interval(Φ * set(F[k - 1]).dat + U.dat)
        _isdisjoint(X, Iₖ) && break
        Δt += δ
        F[k] = ReachSet(Iₖ, Δt)
        k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k - 1)
    end
    return F
end
