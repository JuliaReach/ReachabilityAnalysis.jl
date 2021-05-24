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
                          recursive::Val{false},
                          Δt0::TN) where {N, TN, VNC, VNR}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # preallocations
    n = dim(Ω0)
    c = [VNC(undef, n) for _ in 1:NSTEPS]
    r = [VNR(undef, n) for _ in 1:NSTEPS]

    copy!(c[1], Ω0.center)
    copy!(r[1], Ω0.radius)

    # cache for powers of Φ
    Φ_power_k = similar(Φ)
    copyto!(Φ_power_k, Φ)
    Φ_power_k_cache = similar(Φ_power_k)
    Φ_power_k_abs = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
        # compute the overapproximation of Φ^k * Ω0 with a hyperrectangle
        mul!(c[k], Φ_power_k, c[1])
        Φ_power_k_abs .= abs.(Φ_power_k)
        mul!(r[k], Φ_power_k_abs, r[1])
        Hₖ = Hyperrectangle(c[k], r[k], check_bounds=false)

        Δt += δ
        F[k] = ReachSet(Hₖ, Δt)

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
                          recursive::Val{true},
                          Δt0::TN) where {N, TN, VNC, VNR}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
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
                          recursive::Val{false},
                          Δt0::TN) where {N, TN, VNC, VNR}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # preallocations
    n = dim(Ω0)
    r = [VNR(undef, n) for _ in 1:NSTEPS]
    c = [VNC(undef, n) for _ in 1:NSTEPS]

    copy!(c[1], Ω0.center)
    copy!(r[1], Ω0.radius)

    # cache for powers of Φ
    Φ_power_k = similar(Φ)
    copyto!(Φ_power_k, Φ)
    Φ_power_k_cache = similar(Φ_power_k)
    Φ_power_k_abs = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
      # compute the overapproximation of Φ^k * Ω0 with a hyperrectangle
      mul!(c[k], Φ_power_k, c[1])
      Φ_power_k_abs .= abs.(Φ_power_k)
      mul!(r[k], Φ_power_k_abs, r[1])
      Hₖ = Hyperrectangle(c[k], r[k], check_bounds=false)

      _is_intersection_empty(X, Hₖ) && break
      Δt += δ
      F[k] = ReachSet(Hₖ, Δt)

      mul!(Φ_power_k_cache, Φ_power_k, Φ)
      copyto!(Φ_power_k, Φ_power_k_cache)
      k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k-1)
    end
    return F
end
