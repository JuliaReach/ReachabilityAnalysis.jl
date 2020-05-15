# ===================
# Inhomogeneous case
# ===================

#using LinearAlgebra.BLAS

# inhomogeneous case without invariant, non-recursive implementation
function reach_inhomog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                            Ω0::Hyperrectangle{N, VNC, VNR},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            X::Universe,
                            U::Hyperrectangle,
                            recursive::Val{false},
                            time_shift::N) where {N, VNC, VNR}

    # initial reach set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # preallocations
    c = Vector{VNC}(undef, NSTEPS)
    r = Vector{VNR}(undef, NSTEPS)

    c[1] = Ω0.center
    r[1] = Ω0.radius

    # input sequence
    Wk₊ = copy(U)

    # cache for powers of Φ
    Φ_power_k = _get_cache(Φ)
    Φ_power_k_cache = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
        # Hₖ = overapproximate(Φ_power_k * Ω0, Hyperrectangle)
        abs_Φ_power_k = abs.(Φ_power_k)
        c[k] = Φ_power_k * c[1] + Wk₊.center
        r[k] = abs_Φ_power_k * r[1] + Wk₊.radius
        Hₖ = Hyperrectangle(c[k], r[k], check_bounds=false)

        Δt += δ
        F[k] = ReachSet(Hₖ, Δt)

        # in-place computation of Wk₊ <- Wk₊ + Φ_power_k * U,
        # by overapproximating the second term with a hyperrectangle
        Wk₊.center .+= Φ_power_k * U.center
        Wk₊.radius .+= abs_Φ_power_k * U.radius
        #gemv!('N', 1.0, Φ_power_k, U.center, 1.0, Wk₊.center)
        #gemv!('N', 1.0, abs_Φ_power_k, U.radius, 1.0, Wk₊.radius)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end

# homogeneous case with an invariant, non-recursive implementation
function reach_inhomog_BOX!(F::Vector{ReachSet{N, Hyperrectangle{N, VNC, VNR}}},
                            Ω0::Hyperrectangle{N, VNC, VNR},
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            X::LazySet,
                            U::Hyperrectangle,
                            recursive::Val{false},
                            time_shift::N) where {N, VNC, VNR}

    # initial reach set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = ReachSet(Ω0, Δt)

    # preallocations
    c = Vector{VNC}(undef, NSTEPS)
    r = Vector{VNR}(undef, NSTEPS)

    c[1] = Ω0.center
    r[1] = Ω0.radius

    # input sequence
    Wk₊ = copy(U)

    # cache for powers of Φ
    Φ_power_k = _get_cache(Φ)
    Φ_power_k_cache = similar(Φ_power_k)

    k = 2
    @inbounds while k <= NSTEPS
      # Hₖ = overapproximate(Φ_power_k * Ω0, Hyperrectangle)
      abs_Φ_power_k = abs.(Φ_power_k)
      c[k] = Φ_power_k * c[1] + Wk₊.center
      r[k] = abs.(Φ_power_k) * r[1] + Wk₊.radius
      Hₖ = Hyperrectangle(c[k], r[k], check_bounds=false)

      _is_intersection_empty(X, Hₖ) && break

      Δt += δ
      F[k] = ReachSet(Hₖ, Δt)

      # in-place computation of Wk₊ <- Wk₊ + Φ_power_k * U,
      # by overapproximating the second term with a hyperrectangle
      Wk₊.center .+= Φ_power_k * U.center
      Wk₊.radius .+= abs_Φ_power_k * U.radius

      mul!(Φ_power_k_cache, Φ_power_k, Φ)
      copyto!(Φ_power_k, Φ_power_k_cache)
      k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k-1)
    end
    return F
end
