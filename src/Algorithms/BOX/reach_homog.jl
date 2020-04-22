# assumes that typeof(Φ) is mutable
function _get_cache(Φ::AbstractMatrix)
    return copy(Φ)
end

# change to mutable type
function _get_cache(Φ::SMatrix{S1, S2, T, L}) where {S1, S2, T, L}
    return MMatrix{S1, S2, T, L}(Φ)
end

# ================
# Homogeneous case
# ================

# homogeneous case without invariant, non-recursive implementation
function reach_homog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                          Ω0::Hyperrectangle{N, VNC, VNR},
                          Φ::AbstractMatrix,
                          NSTEPS::Integer,
                          δ::N,
                          X::Universe,
                          recursive::Val{false}) where {N, VNC, VNR}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # cache for powers of Φ
    Φ_power_k = _get_cache(Φ)
    Φ_power_k_cache = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = overapproximate(Φ_power_k * Ω0, Hyperrectangle)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end

# homogeneous case without invariant, recursive implementation
# note that this implementation suffers from wrapping
function reach_homog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                          Ω0::Hyperrectangle{N, VNC, VNR},
                          Φ::AbstractMatrix,
                          NSTEPS::Integer,
                          δ::Float64,
                          X::Universe,
                          recursive::Val{true}) where {N, VNC, VNR}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = overapproximate(Φ * set(F[k-1]), Hyperrectangle)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end

# homogeneous case with an invariant, non-recursive implementation
function reach_homog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                          Ω0::Hyperrectangle{N, VNC, VNR},
                          Φ::AbstractMatrix,
                          NSTEPS::Integer,
                          δ::N,
                          X::LazySet,
                          recursive::Val{false}) where {N, VNC, VNR}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # cache for powers of Φ
    Φ_power_k = _get_cache(Φ)
    Φ_power_k_cache = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = overapproximate(Φ_power_k * Ω0, Hyperrectangle)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end

# homogeneous case with an invariant, recursive implementation
function reach_homog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                          Ω0::Hyperrectangle{N, VNC, VNR},
                          Φ::AbstractMatrix,
                          NSTEPS::Integer,
                          δ::N,
                          X::LazySet,
                          recursive::Val{true}) where {N, VNC, VNR}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = overapproximate(Φ * set(F[k-1]), Hyperrectangle)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end
