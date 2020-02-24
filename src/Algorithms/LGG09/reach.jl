function reach_dir(ρvec_ℓ, N::Int, ℓ::AbstractVector{T},
                   Φₜ::AbstractMatrix{T}, Ω₀::LazySet{T}, Wτ::LazySet{T}) where {T}
    # preallocate output sequence
    ρvec_ℓ = Vector{T}(undef, N-1)
    ρvec = fill(ρvec_ℓ, ndirs)
    @inbounds for (i, ℓ) in enumerate(dirs)
        # note that each of these computations is independent
        reach_dir!(ρvec[i], N, ℓ, Φₜ, Ω₀, Wτ)
    end
end

# TODO: specialize for static vector  static matrix Φₜ ?

# compute N steps of the support function along direction ℓ
function reach_dir!(ρvec_ℓ, N::Int, ℓ::AbstractVector{T},
                   Φₜ::AbstractMatrix{T}, Ω₀::LazySet{T}, Wτ::LazySet{T}) where {T}
    Φₜᵀ = copy(transpose(Φₜ))
    #Φₜᵀ = transpose(Φₜ)
    sᵢ = zero(T) # initialize s₀
    rᵢ = copy(ℓ) # initialize r₀
    rᵢ₊₁ = similar(rᵢ) # initialize r₁
    @inbounds for i in 0:N-2
        mul!(rᵢ₊₁, Φₜᵀ, rᵢ)
        sᵢ₊₁ = sᵢ + ρ(rᵢ, Wτ)
        ρvec_ℓ[i+1] = ρ(rᵢ₊₁, Ω₀) + sᵢ₊₁

        #copy!(rᵢ, rᵢ₊₁) # update cache for the next iteration
        rᵢ .= rᵢ₊₁
        sᵢ = sᵢ₊₁
    end
    return ρvec_ℓ
end
