# ================
# Homogeneous case
# ================

# homogeneous case; no invariant; recursive implementation
function reach_homog_INT!(F::Vector{ReachSet{N, Interval{N, IA.Interval{N}}}},
                          Ω0::Interval{N, IA.Interval{N}},
                          Φ::N,
                          NSTEPS::Integer,
                          δ::Float64,
                          X::Universe) where {N}
    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Iₖ = Interval(Φ * set(F[k-1]).dat)
        Δt += δ
        F[k] = ReachSet(Iₖ, Δt)
        k += 1
    end
    return F
end

# homogeneous case; with invariant; recursive implementation
function reach_homog_INT!(F::Vector{ReachSet{N, Interval{N, IA.Interval{N}}}},
                          Ω0::Interval{N, IA.Interval{N}},
                          Φ::N,
                          NSTEPS::Integer,
                          δ::Float64,
                          X::LazySet) where {N}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Iₖ = Interval(Φ * set(F[k-1]).dat)
        _is_intersection_empty(X, Iₖ) && break
        Δt += δ
        F[k] = ReachSet(Iₖ, Δt)
        k += 1
    end
    if k < NSTEPS +1
        resize!(F, k-1)
    end
    return F
end
