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
                              Δt0::TimeInterval,
                              cache,
                              threaded) where {N,VN,TN,SN,RT<:TemplateReachSet{N,VN,TN,SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ρℓ = Matrix{N}(undef, length(dirs), NSTEPS)

    _reach_inhomog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, U, dirs, NSTEPS, cache, threaded)

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + Δt0
    @inbounds for k in 1:NSTEPS
        F[k] = TemplateReachSet(dirs, view(ρℓ, :, k), Δt)
        Δt += δ
    end

    return ρℓ
end

function _reach_inhomog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, U, dirs, NSTEPS, cache, threaded::Val{true})
    ℓ = _collect(dirs)
    Threads.@threads for j in eachindex(ℓ)
        reach_inhomog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, U, ℓ[j], NSTEPS, cache)
    end
end

function _reach_inhomog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, U, dirs, NSTEPS, cache, threaded::Val{false})
    #for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, U, ℓ, NSTEPS, cache)
    end
end

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, U, ℓ::AbstractVector{N},
                                  NSTEPS, cache::Val{true}) where {N}
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

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, U, ℓ::AbstractVector{N},
                                  NSTEPS, cache::Val{false}) where {N}
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

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, U, ℓ::AbstractVector{N},
                                  NSTEPS, cache::Val{true}) where {N}
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

function reach_inhomog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, U, ℓ::AbstractVector{N},
                                  NSTEPS, cache::Val{false}) where {N}
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
    ρℓ = [Vector{N}(undef, NSTEPS) for _ in dirs]

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_inhomog_dir_LGG09!(ρℓ[i], Ω₀, Φᵀ, U, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        sf = [ρℓ[i][k] for i in eachindex(dirs)]
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

# this version uses one matrix-matrix product for each time-step
function reach_inhomog_LGG09!(F::Vector{RT},
                              dirs::TN,
                              Ω₀::LazySet{N},
                              Φ::AbstractMatrix{N},
                              NSTEPS::Integer,
                              δ::N,
                              X::LazySet,
                              U::LazySet,
                              Δt0::TimeInterval,
                              cache,
                              threaded) where {N,VN,TN,SN,RT<:TemplateReachSet{N,VN,TN,SN}}

    # transpose coefficients matrix
    Φᵀ = copy(transpose(Φ))

    # preallocate output sequence
    ndirs = length(dirs)
    ρmat = Matrix{N}(undef, ndirs, NSTEPS)

    rᵢ = reduce(hcat, dirs) # one column per template direction
    rᵢ₊₁ = similar(rᵢ)
    sᵢ = zeros(N, ndirs)

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + Δt0
    k = 1
    @inbounds while k <= NSTEPS

        # loop dispatch in threaded
        _reach_inhomog_LGG09_invariant!(Ω₀, U, rᵢ, ρmat, sᵢ, k, ndirs, threaded)

        # update cache for the next iteration
        mul!(rᵢ₊₁, Φᵀ, rᵢ)
        copy!(rᵢ, rᵢ₊₁)

        F[k] = TemplateReachSet(dirs, view(ρmat, :, k), Δt)
        _is_intersection_empty(X, set(F[k])) && break  # TODO pass disjointness method
        Δt += δ
        k += 1
    end
    if k < NSTEPS
        resize!(F, k - 1)
    end
    return ρmat
end

function _reach_inhomog_LGG09_invariant!(Ω₀, U, rᵢ, ρmat, sᵢ, k, ndirs, threaded::Val{false})
    for j in 1:ndirs
        d = view(rᵢ, :, j)
        ρmat[j, k] = ρ(d, Ω₀) + sᵢ[j]
        sᵢ[j] += ρ(d, U)
    end
end

function _reach_inhomog_LGG09_invariant!(Ω₀, U, rᵢ, ρmat, sᵢ, k, ndirs, threaded::Val{true})
    Threads.@threads for j in 1:ndirs
        d = view(rᵢ, :, j)
        ρmat[j, k] = ρ(d, Ω₀) + sᵢ[j]
        sᵢ[j] += ρ(d, U)
    end
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

# ------------------------------------------------------------
# Functionality that requires ExponentialUtilities.jl
# ------------------------------------------------------------

function load_krylov_LGG09_inhomog()
    return quote

        # Compute the sequence with constant input sets:
        #
        # ρ(ℓ, Ω₀), ρ(ℓ, Φ Ω₀ ⊕ V), ρ(ℓ, Φ^2 Ω₀ ⊕ Φ V ⊕ V), ρ(ℓ, Φ^3 Ω₀ ⊕ Φ^2 V ⊕ Φ V ⊕ V), ...
        #
        # Using Krylov subspace approximations to compute the action of Φ := exp(Aδ) over
        # the direction ℓ.
        #
        # Method (see [1]):
        #
        # out[1] <- ρ(ℓ, Ω₀)
        #
        # out[2] <- ρ(ℓ, Φ Ω₀ ⊕ V) = ρ(ℓ, Φ Ω0) + ρ(ℓ, V) = ρ(Φᵀ ℓ, Ω₀) + ρ(ℓ, V)
        #
        # out[3] <- ρ(ℓ, Φ^2 Ω₀ ⊕ Φ V ⊕ V) = ρ((Φᵀ)^2 ℓ, Ω₀) + ρ(Φᵀ ℓ, V) + ρ(ℓ, V)
        #
        # out[4] <- ρ(ℓ, Φ^3 Ω₀ ⊕ Φ^2 V ⊕ Φ V ⊕ V) = ρ((Φᵀ)^3 ℓ, Ω₀) + ρ((Φᵀ)^2 ℓ, V) + ρ(Φᵀ ℓ, V) + ρ(ℓ, V)
        #
        # and so on.
        #
        # [1] Reach Set Approximation through Decomposition with Low-dimensional Sets and
        #     High-dimensional Matrices. Sergiy Bogomolov, Marcelo Forets, Goran Frehse,
        #     Frédéric Viry, Andreas Podelski and Christian Schilling (2018) HSCC'18
        #     Proceedings of the 21st International Conference on Hybrid Systems: Computation
        #     and Control: 41–50.
        #
        function reach_inhomog_krylov_LGG09!(out, Ω₀::LazySet, V::LazySet, Aᵀδ::AbstractMatrix,
                                             ℓ::AbstractVector, NSTEPS;
                                             hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

            # initialization of the krylov subspace
            TA, Tb = eltype(Aᵀδ), eltype(ℓ)
            T = promote_type(TA, Tb)
            Ks = KrylovSubspace{T,real(T)}(length(ℓ), m)
            arnoldi!(Ks, Aᵀδ, ℓ; m=m, ishermitian=hermitian, tol=tol)

            # rᵢ stores is the cache for each vector: (Φᵀ)^i ℓ
            rᵢ = deepcopy(ℓ)
            out[1] = ρ(ℓ, Ω₀)

            # accumulated support vector sum due to the inputs
            s = zero(T)

            @inbounds for i in 1:(NSTEPS - 1)
                s += ρ(rᵢ, V)
                expv!(rᵢ, i * 1.0, Ks)
                out[i + 1] = ρ(rᵢ, Ω₀) + s
            end
            return out
        end

        # Compute the sequence with time-varying input sets:
        #
        # ρ(ℓ, Ω₀), ρ(ℓ, Φ Ω₀ ⊕ V₀), ρ(ℓ, Φ^2 Ω₀ ⊕ Φ V₀ ⊕ V₁), ρ(ℓ, Φ^3 Ω₀ ⊕ Φ^2 V₀ ⊕ Φ V₁ ⊕ V₂), ...
        #
        # Using Krylov subspace approximations to compute the action of Φ := exp(Aδ) over
        # the direction ℓ.
        #
        # Method (see [1]):
        #
        # out[1] <- ρ(ℓ, Ω₀)
        #
        # out[2] <- ρ(ℓ, Φ Ω₀ ⊕ V₀) = ρ(ℓ, Φ Ω0) + ρ(ℓ, V₀) = ρ(Φᵀ ℓ, Ω₀) + ρ(ℓ, V₀)
        #
        # out[3] <- ρ(ℓ, Φ^2 Ω₀ ⊕ Φ V₀ ⊕ V₁) = ρ((Φᵀ)^2 ℓ, Ω₀) + ρ(Φᵀ ℓ, V₀) + ρ(ℓ, V₁)
        #
        # out[4] <- ρ(ℓ, Φ^3 Ω₀ ⊕ Φ^2 V₀ ⊕ Φ V₁ ⊕ V₂) = ρ((Φᵀ)^3 ℓ, Ω₀) + ρ((Φᵀ)^2 ℓ, V₀) + ρ(Φᵀ ℓ, V₁) + ρ(ℓ, V₂)
        #
        # and so on.
        #
        # [1] Reach Set Approximation through Decomposition with Low-dimensional Sets and
        #     High-dimensional Matrices. Sergiy Bogomolov, Marcelo Forets, Goran Frehse,
        #     Frédéric Viry, Andreas Podelski and Christian Schilling (2018) HSCC'18
        #     Proceedings of the 21st International Conference on Hybrid Systems: Computation
        #     and Control: 41–50.
        #
        function reach_inhomog_krylov_LGG09!(out, Ω₀::LazySet, V::Vector{<:LazySet},
                                             Aᵀδ::AbstractMatrix,
                                             ℓ::AbstractVector, NSTEPS;
                                             hermitian=false, m=min(30, size(Aᵀδ, 1)), tol=1e-7)

            # initialization of the krylov subspace
            TA, Tb = eltype(Aᵀδ), eltype(ℓ)
            T = promote_type(TA, Tb)
            Ks = KrylovSubspace{T,real(T)}(length(ℓ), m)
            arnoldi!(Ks, Aᵀδ, ℓ; m=m, ishermitian=hermitian, tol=tol)

            # rᵢ stores is the cache for each vector: (Φᵀ)^i ℓ
            r = [similar(ℓ) for _ in NSTEPS]
            r[1] = deepcopy(ℓ)
            out[1] = ρ(ℓ, Ω₀)

            @inbounds for i in 1:(NSTEPS - 1)
                s = zero(T)
                for j in 1:i
                    s += ρ(r[i - j + 1], V[j])
                end
                expv!(r[i + 1], i * 1.0, Ks)
                out[i + 1] = ρ(r[i + 1], Ω₀) + s
            end
            return out
        end
    end
end  # quote / load_krylov_LGG09_inhomog()

# ------------------------------------------------------------
# Methods using eigenvalues of the transition matrix
# ------------------------------------------------------------

#= TODO needs review <<

# it is assumed that (λ, d) is an eigenvalue-eigenvector pair of the matrix Φᵀ
function reach_inhomog_dir_eig_LGG09!(out::AbstractVector{N}, X₀, V, d::AbstractVector{N}, λ::N, NSTEPS) where {N}
    if iszero(λ)
        _reach_inhomog_dir_eig_LGG09_zero!(out, V, NSTEPS)
    elseif λ > zero(N)
        _reach_inhomog_dir_eig_LGG09_positive!(out, X₀, V, d, λ, NSTEPS)
    else
        _reach_inhomog_dir_eig_LGG09_negative!(out, X₀, V, d, λ, NSTEPS)
    end
    return out
end

function _reach_inhomog_dir_eig_LGG09_zero!(out::AbstractVector{N}, V::LazySet, NSTEPS) where {N}
    ρ_d_V = ρ(d, V)
    @inbounds for i in 1:NSTEPS
        out[i] = ρ_d_V
    end
end

function _reach_inhomog_dir_eig_LGG09_zero!(out::AbstractVector{N}, V::AbstractVector{<:LazySet}, NSTEPS) where {N}
    @inbounds for i in 1:NSTEPS
        out[i] = ρ(d, V[i])
    end
end

function _reach_inhomog_dir_eig_LGG09_positive!(out::AbstractVector{N}, X₀, V::LazySet, d, λ, NSTEPS) where {N}
    ρ_d_X₀ = ρ(d, X₀)
    ρ_d_V = ρ(d, V)
    @inbounds begin
        out[1] = ρ_d_X₀
        λⁱ = λ
        λacc = one(N)
        for i in 2:NSTEPS
            out[i] = λⁱ * ρ_d_X₀ + λacc * ρ_d_V
            λacc = λacc + λⁱ
            λⁱ = λⁱ * λ
        end
    end
end

function _reach_inhomog_dir_eig_LGG09_positive!(out::AbstractVector{N}, X₀, V::AbstractVector{<:LazySet}, d, λ, NSTEPS) where {N}
    ρ_d_X₀ = ρ(d, X₀)
    @inbounds begin
        out[1] = ρ_d_X₀
        ρ_d_V_acc = ρ(d, V[1])
        λⁱ = λ
        for i in 2:NSTEPS
            out[i] = λⁱ * ρ_d_X₀ + ρ_d_V_acc
            ρ_d_V_acc = λ * ρ_d_V_acc + ρ(d, V[i])
            λⁱ = λⁱ * λ
        end
    end
end

function _reach_inhomog_dir_eig_LGG09_negative!(out::AbstractVector{N}, X₀, V::LazySet, d, λ, NSTEPS) where {N}
    ρ_d_X₀ = ρ(d, X₀)
    ρ_md_X₀ = ρ(-d, X₀)
    ρ_d_V = ρ(d, V)
    ρ_md_V = ρ(-d, V)
    λⁱ = λ
    @inbounds begin
        out[1] = ρ_d_X₀
        for i in 2:NSTEPS
            if iseven(i)
                out[i] = -λⁱ * ρ_md_X₀
            else
                out[i] = λⁱ * ρ_d_X₀
            end
            λⁱ = λⁱ * λ
        end
    end
end

function _reach_inhomog_dir_eig_LGG09_negative!(out::AbstractVector{N}, X₀, V::AbstractVector{<:LazySet}, d, λ, NSTEPS) where {N}
    ρ_d_X₀ = ρ(d, X₀)
    ρ_md_X₀ = ρ(-d, X₀)
    λⁱ = λ
    @inbounds begin
        out[1] = ρ_d_X₀
        for i in 2:NSTEPS
            if iseven(i)
                out[i] = -λⁱ * ρ_md_X₀
            else
                out[i] = λⁱ * ρ_d_X₀
            end
            λⁱ = λⁱ * λ
        end
    end
end
=#
