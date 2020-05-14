# ===================
# Inhomogeneous case
# ===================

# note that each of the computations in the loop over directions is independent
function reach_inhomog_LGG09!(ρvec_ℓ, Φ::AbstractMatrix{N}, Ω₀::LazySet{N},
                              ℓ::AbstractVector{N}, NSTEPS::Int, W::LazySet{N}, t0) where {N}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρvec_ℓ = Vector{T}(undef, N-1)
    ρvec = fill(ρvec_ℓ, ndirs)

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρvec[i], Φᵀ, ℓ, Ω₀, NSTEPS, W)
    end
end

# TODO: needs specialization for static vector / static matrix ?

# compute NSTEPS iterations support function along direction ℓ
function reach_inhomog_dir_LGG09!(ρvec_ℓ, Φᵀ, ℓ, Ω₀, NSTEPS, W, t0)
    sᵢ = zero(T) # initialize s₀
    rᵢ = copy(ℓ) # initialize r₀
    rᵢ₊₁ = similar(rᵢ) # initialize r₁

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
