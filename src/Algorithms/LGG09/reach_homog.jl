# =================
# Homogeneous case
# =================

# note that each of the computations in the loop over directions is independent
function reach_homog_LGG09!(F::Vector{TemplateReachSet{N, VN, TN, SN}}, # SN = Vector{N} always ?
                            dirs::TN,
                            Ω₀::LazySet,
                            Φ::AbstractMatrix,
                            NSTEPS::Integer,
                            δ::N,
                            X::Universe, # no invariant
                            time_shift::N
                            ) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    # for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return F
end

function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, ℓ, NSTEPS) where {N}
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

# TODO: needs specialization for static vector / static matrix ?
# compute NSTEPS iterations support function along direction ℓ
function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, ℓ, NSTEPS) where {N}
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
