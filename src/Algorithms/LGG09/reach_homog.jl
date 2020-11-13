# =================
# Homogeneous case
# =================

# note that each of the computations in the loop over directions is independent
function reach_homog_LGG09!(F::Vector{RT},
                            dirs::TN,
                            Ω₀::LazySet,
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            X::Universe, # no invariant
                            Δt0::TimeInterval,
                            cache,
                            threaded) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    _reach_homog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, dirs, NSTEPS, cache, threaded)

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + Δt0
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return ρℓ
end

function _reach_homog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, dirs, NSTEPS, cache, threaded::Val{true})
    ℓ = _collect(dirs)
    Threads.@threads for j in 1:length(ℓ)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ[j], NSTEPS, cache)
    end
end

function _reach_homog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, dirs, NSTEPS, cache, threaded::Val{false})
    #for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ, NSTEPS, cache)
    end
end

function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, ℓ, NSTEPS, cache::Val{true}) where {N}
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[j, i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)   # rᵢ .= rᵢ₊₁
    end
    return ρvec_ℓ
end

function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, ℓ, NSTEPS, cache::Val{false}) where {N}
    rᵢ = copy(ℓ)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[j, i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        rᵢ = Φᵀ * rᵢ
    end
    return ρvec_ℓ
end

# TODO: needs specialization for static vector / static matrix ?
# compute NSTEPS iterations support function along direction ℓ
# ``ρ(ℓ, Ω₀), ρ(ℓ, Φ * Ω₀), ρ(ℓ, Φ^2 * Ω₀), ..., ρ(ℓ, Φ^NSTEPS * Ω₀)``
function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, ℓ, NSTEPS, cache::Val{true}) where {N}
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)   # rᵢ .= rᵢ₊₁
    end
    return ρvec_ℓ
end

# ``ρ(ℓ, Ω₀), ρ(ℓ, Φ * Ω₀), ρ(ℓ, Φ^2 * Ω₀), ..., ρ(ℓ, Φ^NSTEPS * Ω₀)``
function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, ℓ, NSTEPS, cache::Val{false}) where {N}
    rᵢ = copy(ℓ)
    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        rᵢ = Φᵀ * rᵢ
    end
    return ρvec_ℓ
end

# NOTES:
# - if ℓ and -ℓ are in the template, they could share (Φᵀ)^k * ℓ
# - the loop over directions can be parallelized

#=---- version using vector-of-vectors
    # preallocate output sequence
    ρℓ = [Vector{N}(undef, NSTEPS) for _ in 1:length(dirs)]

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ[i], Ω₀, Φᵀ, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        sf = [ρℓ[i][k] for i in 1:length(dirs)]
        F[k] = TemplateReachSet(dirs, sf, Δt)
        Δt += δ
    end
=#

# ------------------------------------------------------------
# Functionality that requires ExponentialUtilities.jl
# ------------------------------------------------------------

function load_krylov_LGG09_homog()
return quote

"""
    reach_homog_krylov_LGG09!(out, Ω₀::LazySet, Aᵀδ::AbstractMatrix,
                              ℓ::AbstractVector, NSTEPS;
                              hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

### Algorithm

We compute the sequence:

```math
    ρ(ℓ, Ω₀), ρ(ℓ, Φ Ω₀), ρ(ℓ, Φ^2 Ω₀), ρ(ℓ, Φ^3 Ω₀), ...
```

Using Krylov subspace approximations to compute the action of Φ := exp(Aδ) over
the direction ℓ.

The method is (see [1]):

```julia
out[1] <- ρ(ℓ, Ω₀)

out[2] <- ρ(ℓ, Φ Ω₀) = ρ(Φᵀ ℓ, Ω₀)

out[3] <- ρ(ℓ, Φ^2 Ω₀) = ρ((Φᵀ)^2 ℓ, Ω₀)

out[4] <- ρ(ℓ, Φ^3 Ω₀) = ρ((Φᵀ)^3 ℓ, Ω₀)
```
 and so on.

### References

[1] Reach Set Approximation through Decomposition with Low-dimensional Sets and
    High-dimensional Matrices. Sergiy Bogomolov, Marcelo Forets, Goran Frehse,
    Frédéric Viry, Andreas Podelski and Christian Schilling (2018) HSCC'18
    Proceedings of the 21st International Conference on Hybrid Systems: Computation
    and Control: 41–50.
"""
function reach_homog_krylov_LGG09!(out, Ω₀::LazySet, Aᵀδ::AbstractMatrix,
                                   ℓ::AbstractVector, NSTEPS;
                                   hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

    # initialization of the krylov subspace
    TA, Tb = eltype(Aᵀδ), eltype(ℓ)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(ℓ), m)
    arnoldi!(Ks, Aᵀδ, ℓ; m=m, ishermitian=hermitian, tol=tol)

    # rᵢ stores is the cache for each vector: (Φᵀ)^i ℓ
    rᵢ = deepcopy(ℓ)

    @inbounds for i in 1:NSTEPS
        out[i] = ρ(rᵢ, Ω₀)
        expv!(rᵢ, i*1.0, Ks)
    end
    return out
end

end end  # quote / load_krylov_LGG09_homog()
