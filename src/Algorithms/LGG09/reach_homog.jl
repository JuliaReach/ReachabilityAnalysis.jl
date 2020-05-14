# =================
# Homogeneous case
# =================

# note that each of the computations in the loop over directions is independent
function reach_homog_LGG09!(ρvec_ℓ, Φ::AbstractMatrix{N}, Ω₀::LazySet{N},
                            ℓ::AbstractVector{N}, NSTEPS::Int, t0) where {N}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρvec_ℓ = Vector{N}(undef, N-1)
    ρvec = fill(ρvec_ℓ, ndirs)

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρvec[i], Φᵀ, ℓ, Ω₀, NSTEPS)
    end
    return ρvec_ℓ
end

# TODO: needs specialization for static vector / static matrix ?

# compute NSTEPS iterations support function along direction ℓ
function reach_homog_dir_LGG09!(ρvec_ℓ, Φᵀ, ℓ, Ω₀, NSTEPS, t0) where {N}
    rᵢ = copy(ℓ) # initialize r₀
    rᵢ₊₁ = similar(rᵢ) # initialize r₁

    @inbounds for i in 0:NSTEPS-2
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        ρvec_ℓ[i+1] = ρ(rᵢ₊₁, Ω₀)

        #copy!(rᵢ, rᵢ₊₁) # update cache for the next iteration
        rᵢ .= rᵢ₊₁
    end
    return ρvec_ℓ
end
