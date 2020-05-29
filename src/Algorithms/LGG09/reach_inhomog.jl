# ===================
# Inhomogeneous case
# ===================

function reach_inhomog_LGG09!(F::Vector{RT},
                              dirs::TN,
                              Ω₀::LazySet{N},
                              Φ::AbstractMatrix{N},
                              NSTEPS::Integer,
                              δ::N,
                              X::Universe,
                              U::LazySet,
                              time_shift::N) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    # for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, U, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return ρℓ
end

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS) where {N}
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)
    sᵢ = zero(N)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[j, i] = ρ(rᵢ, Ω₀) + sᵢ
        sᵢ += ρ(rᵢ, U)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)   # rᵢ .= rᵢ₊₁
    end
    return ρvec_ℓ
end

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS) where {N}
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)
    sᵢ = zero(N)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[i] = ρ(rᵢ, Ω₀) + sᵢ
        sᵢ += ρ(rᵢ, U)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)   # rᵢ .= rᵢ₊₁
    end
    return ρvec_ℓ
end

#= ---- version using vector-of-vectors
    ρℓ = [Vector{N}(undef, NSTEPS) for _ in 1:length(dirs)]

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ[i], Ω₀, Φᵀ, U, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        sf = [ρℓ[i][k] for i in 1:length(dirs)]
        F[k] = TemplateReachSet(dirs, sf, Δt)
        Δt += δ
    end
---- =#

#=
# TODO: specialization for static vector / static matrix

# compute NSTEPS iterations support function along direction ℓ
function reach_inhomog_dir_LGG09!(ρvec_ℓ, Φᵀ, ℓ, Ω₀, NSTEPS, W, time_shift)
    sᵢ = zero(T)
    rᵢ = copy(ℓ)
    rᵢ₊₁ = similar(rᵢ)

    @inbounds for i in 0:NSTEPS-2
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        sᵢ₊₁ = sᵢ + ρ(rᵢ, W)
        ρvec_ℓ[i+1] = ρ(rᵢ₊₁, Ω₀) + sᵢ₊₁

        #copy!(rᵢ, rᵢ₊₁) # update cache for the next iteration
        rᵢ .= rᵢ₊₁
        sᵢ = sᵢ₊₁
    end
    return ρvec_ℓ
end
=#
