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
                            time_shift::N,
                            cache) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    # for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ, NSTEPS, cache)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return ρℓ
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

function load_exponential_utilities_LGG09()
return quote

# recursive version, default expv
function reach_homog_dir_LGG09_expv!(out, Ω₀, Aᵀ, ℓ, NSTEPS, recursive::Val{:true})
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 1:NSTEPS
        out[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        rᵢ₊₁ = expv(1.0, Aᵀ, rᵢ) # computes exp(Aᵀ * 1.0) * rᵢ
        copy!(rᵢ, rᵢ₊₁)
    end
    return out
end

# ρ(ℓ, Ω₀), ρ(exp(Aᵀ) * ℓ, Ω₀), ρ(exp(2Aᵀ) * ℓ, Ω₀)
function reach_homog_dir_LGG09_expv!(out, Ω₀, Aᵀ, ℓ, NSTEPS, recursive::Val{:false})
    rᵢ = deepcopy(ℓ) # if ℓ is a sev => this is a sev

    @inbounds for i in 1:NSTEPS
        out[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        rᵢ = expv(i*1.0, Aᵀ, ℓ)
    end
    return out
end

# non-recursive version using precomputed Krylov subspace
function reach_homog_dir_LGG09_expv_pk!(out, Ω₀, Aᵀ, ℓ, NSTEPS, recursive::Val{:false})
    rᵢ = deepcopy(ℓ) # if ℓ is a sev => this is a sev
    Ks = arnoldi(Aᵀ, ℓ, tol=1e-18)

    @inbounds for i in 1:NSTEPS
        out[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        expv!(rᵢ, i*1.0, Ks)
    end
    return out
end

# this function computes the sequence
# ``ρ(ℓ, Ω₀)``, ``ρ(exp(Aᵀ) * ℓ, Ω₀)``, ``ρ(exp(2Aᵀ) * ℓ, Ω₀)`` until ``ρ(exp(NSTEPS * Aᵀ) * ℓ, Ω₀)``
function reach_homog_dir_LGG09_expv_pk2!(out, Ω₀, Aᵀ, ℓ, NSTEPS, recursive::Val{:false};
                                        hermitian=false, m=min(30, size(Aᵀ, 1)), tol=1e-7)

    TA, Tb = eltype(Aᵀ), eltype(ℓ)
    T = promote_type(TA, Tb)
    Ks = KrylovSubspace{T, real(T)}(length(ℓ), m)
    arnoldi!(Ks, Aᵀ, ℓ; m=m, ishermitian=hermitian, tol=tol)

    rᵢ = deepcopy(ℓ)

    @inbounds for i in 1:NSTEPS
        out[i] = ρ(rᵢ, Ω₀)

        # update cache for the next iteration
        expv!(rᵢ, i*1.0, Ks)
    end
    return out
end

end end  # quote / load_exponential_utilities_LGG09()
