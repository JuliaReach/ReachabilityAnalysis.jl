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
                              time_shift::N,
                              cache) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    # for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, U, ℓ, NSTEPS, cache)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return ρℓ
end

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS, cache::Val{true}) where {N}
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

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS, cache::Val{false}) where {N}
    rᵢ = copy(ℓ)
    sᵢ = zero(N)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[j, i] = ρ(rᵢ, Ω₀) + sᵢ
        sᵢ += ρ(rᵢ, U)

        # update cache for the next iteration
        rᵢ = Φᵀ * rᵢ
    end
    return ρvec_ℓ
end

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS, cache::Val{true}) where {N}
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

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, U, ℓ::AbstractVector{N}, NSTEPS, cache::Val{false}) where {N}
    rᵢ = copy(ℓ)
    sᵢ = zero(N)

    @inbounds for i in 1:NSTEPS
        ρvec_ℓ[i] = ρ(rᵢ, Ω₀) + sᵢ
        sᵢ += ρ(rᵢ, U)

        # update cache for the next iteration
        rᵢ = Φᵀ * rᵢ
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


# ===================================
# Inhomogeneous case with invariant
# ===================================

# in this version we one matrix-matrix product for each time-step
function reach_inhomog_LGG09!(F::Vector{RT},
                              dirs::TN,
                              Ω₀::LazySet{N},
                              Φ::AbstractMatrix{N},
                              NSTEPS::Integer,
                              δ::N,
                              X::LazySet,
                              U::LazySet,
                              time_shift::N,
                              cache) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ndirs = length(dirs)
    ρmat = Matrix{N}(undef, ndirs, NSTEPS)

    rᵢ = reduce(hcat, dirs) # one column per template direction
    rᵢ₊₁ = similar(rᵢ)
    sᵢ = zeros(N, ndirs)

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    k = 1
    @inbounds while k <= NSTEPS
        for j in 1:ndirs
            d = view(rᵢ, :, j)
            ρmat[j, k] = ρ(d, Ω₀) + sᵢ[j]
            sᵢ[j] += ρ(d, U)
        end
        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)

        F[k] = TemplateReachSet(dirs, view(ρmat, :, k), Δt)
        _is_intersection_empty(X, set(F[k])) && break  # TODO pass disjointness method
        Δt += δ
        k += 1
    end
    if k < NSTEPS
        resize!(F, k-1)
    end
    return ρmat
end


#=
# in this version we use several matrix-vector products for each direction
function reach_inhomog_LGG09!(F::Vector{RT},
                              dirs::TN,
                              Ω₀::LazySet{N},
                              Φ::AbstractMatrix{N},
                              NSTEPS::Integer,
                              δ::N,
                              X::LazySet,
                              U::LazySet,
                              time_shift::N,
                              cache) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ndirs = length(dirs)
    ρmat = Matrix{N}(undef, ndirs, NSTEPS)

    rᵢ = [copy(ℓ) for ℓ in dirs]
    rᵢ₊₁ = [similar(first(rᵢ)) for _ in 1:ndirs]
    sᵢ = zeros(N, ndirs)

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    k = 1
    @inbounds while k <= NSTEPS
        for (j, ℓ) in enumerate(dirs)
            ρmat[j, k] = ρ(rᵢ[j], Ω₀) + sᵢ[j]
            sᵢ[j] += ρ(rᵢ[j], U)

            # update cache for the next iteration
            mul!(rᵢ₊₁[j], Φᵀ, rᵢ[j])
            copy!(rᵢ[j], rᵢ₊₁[j])
        end

        F[k] = TemplateReachSet(dirs, view(ρmat, :, k), Δt)
        _is_intersection_empty(X, set(F[k])) && break  # TODO pass disjointness method
        Δt += δ
        k += 1
    end
    if k < NSTEPS
        resize!(F, k-1)
    end
    return ρmat
end
=#

#=
# in this version all the reach-sets are computed until NSTEPS, then we intersect
# with the invariants afterwards
function reach_inhomog_LGG09!(F::Vector{RT},
                              dirs::TN,
                              Ω₀::LazySet{N},
                              Φ::AbstractMatrix{N},
                              NSTEPS::Integer,
                              δ::N,
                              X::LazySet,
                              U::LazySet,
                              time_shift::N,
                              cache) where {N, VN, TN, SN, RT<:TemplateReachSet{N, VN, TN, SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    # for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, U, ℓ, NSTEPS, cache)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    k = 1
    @inbounds while k <= NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        _is_intersection_empty(X, set(F[k])) && break  # TODO pass disjointness method
        Δt += δ
        k += 1
    end
    if k < NSTEPS
        resize!(F, k-1)
    end
    return ρℓ
end
=#
