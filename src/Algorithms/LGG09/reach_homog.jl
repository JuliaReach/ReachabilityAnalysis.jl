# ===================================
# Homogeneous case without invariant
# ===================================

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
                            threaded) where {N,VN,TN,SN,RT<:TemplateReachSet{N,VN,TN,SN}}

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
    Threads.@threads for j in eachindex(ℓ)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ[j], NSTEPS, cache)
    end
end

function _reach_homog_dir_LGG09!(ρℓ, Ω₀, Φᵀ, dirs, NSTEPS, cache, threaded::Val{false})
    #for each direction, compute NSTEPS iterations
    @inbounds for (j, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ, j, Ω₀, Φᵀ, ℓ, NSTEPS, cache)
    end
end

function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, ℓ, NSTEPS,
                                cache::Val{true}) where {N}
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

function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractMatrix{N}, j, Ω₀, Φᵀ, ℓ, NSTEPS,
                                cache::Val{false}) where {N}
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
function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, ℓ, NSTEPS,
                                cache::Val{true}) where {N}
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
function reach_homog_dir_LGG09!(ρvec_ℓ::AbstractVector{N}, Ω₀, Φᵀ, ℓ, NSTEPS,
                                cache::Val{false}) where {N}
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
    ρℓ = [Vector{N}(undef, NSTEPS) for _ in dirs]

    # for each direction, compute NSTEPS iterations
    @inbounds for (i, ℓ) in enumerate(dirs)
        reach_homog_dir_LGG09!(ρℓ[i], Ω₀, Φᵀ, ℓ, NSTEPS)
    end

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + time_shift
    @inbounds for k in 1:NSTEPS
        sf = [ρℓ[i][k] for i in eachindex(dirs)]
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
            Ks = KrylovSubspace{T,real(T)}(length(ℓ), m)
            arnoldi!(Ks, Aᵀδ, ℓ; m=m, ishermitian=hermitian, tol=tol)

            # rᵢ stores is the cache for each vector: (Φᵀ)^i ℓ
            rᵢ = deepcopy(ℓ)

            @inbounds for i in 1:NSTEPS
                out[i] = ρ(rᵢ, Ω₀)
                expv!(rᵢ, i * 1.0, Ks)
            end
            return out
        end
    end
end  # quote / load_krylov_LGG09_homog()

# ------------------------------------------------------------
# Methods using eigenvalues of the transition matrix
# ------------------------------------------------------------

# it is assumed that (λ, d) is an eigenvalue-eigenvector pair of the matrix Φᵀ
function reach_homog_dir_eig_LGG09!(out::AbstractVector{N}, X₀, d::AbstractVector{N}, λ::N,
                                    NSTEPS) where {N}
    if iszero(λ)
        _reach_homog_dir_eig_LGG09_zero!(out, X₀, d, NSTEPS)
    elseif λ > zero(N)
        _reach_homog_dir_eig_LGG09_positive!(out, X₀, d, λ, NSTEPS)
    else
        _reach_homog_dir_eig_LGG09_negative!(out, X₀, d, λ, NSTEPS)
    end
    return out
end

function _reach_homog_dir_eig_LGG09_zero!(out::AbstractVector{N}, X₀, d, NSTEPS) where {N}
    @inbounds begin
        out[1] = ρ(d, X₀)
        for i in 2:NSTEPS
            out[i] = zero(N)
        end
    end
end

function _reach_homog_dir_eig_LGG09_positive!(out::AbstractVector{N}, X₀, d, λ, NSTEPS) where {N}
    ρ₀ = ρ(d, X₀)
    λⁱ = one(N)
    @inbounds for i in 1:NSTEPS
        out[i] = λⁱ * ρ₀
        λⁱ = λⁱ * λ
    end
end

function _reach_homog_dir_eig_LGG09_negative!(out::AbstractVector{N}, X₀, d, λ, NSTEPS) where {N}
    ρ₀ = ρ(d, X₀)
    ρ₀₋ = ρ(-d, X₀)
    λⁱ = λ
    @inbounds begin
        out[1] = ρ₀
        for i in 2:NSTEPS
            if iseven(i)
                out[i] = -λⁱ * ρ₀₋
            else
                out[i] = λⁱ * ρ₀
            end
            λⁱ = λⁱ * λ
        end
    end
end

# Let A be an `n x n` matrix with real eigenvalues.
# Given the eigendecomposition: Aδ = P * Λ * P^{-1},
# where Λ is a diagonal matrix that contains the eigenvalues of Aδ
# and P is an invertible matrix.
# Let Q = (P^{-1})^T, and let Φ = exp(Aδ). Then it holds that:
# Φ^T = Q * exp(Λ) * Q^{-1}
# Moreover, the eigenvectors of Φ^T are the columns of Q.
#
# This function returns a matrix `ρmat` of size `2n × NSTEPS` such that
# the i-th row of `ρmat` contains all support function evaluations
# ρ(Φ^k dᵢ, Ω₀) for each k = 0, …, NSTEPS-1,
# where dᵢ = vᵢ and vᵢ is the i-th eigenvector of Φ^T if 1 ≤ i ≤ n and for i > n
# dᵢ = -vᵢ is the negative of the i-th eigenvector of Φ^T.
function reach_homog_eig_LGG09_posneg(Λ::Vector{N}, Q, Ω₀, NSTEPS) where {N}
    n = length(Λ)
    @assert n == size(Q, 1) == size(Q, 2)

    ρmat = Matrix{N}(undef, 2 * n, NSTEPS)
    Q₋ = -Q
    for j in 1:n
        λj = exp(Λ[j])

        Qj₊ = view(Q, :, j)
        reach_homog_dir_eig_LGG09!(view(ρmat, j, :), Ω₀, Qj₊, λj, NSTEPS)

        Qj₋ = view(Q₋, :, j)
        reach_homog_dir_eig_LGG09!(view(ρmat, j + n, :), Ω₀, Qj₋, λj, NSTEPS)
    end
    return ρmat
end

# same as reach_homog_eig_LGG09_posneg, but only computes along Q
function reach_homog_eig_LGG09(Λ::Vector{N}, Q, Ω₀, NSTEPS) where {N}
    n = length(Λ)
    @assert n == size(Q, 1) == size(Q, 2)

    ρmat = Matrix{N}(undef, n, NSTEPS)
    for j in 1:n
        λj = exp(Λ[j])
        Qj = view(Q, :, j)
        reach_homog_dir_eig_LGG09!(view(ρmat, j, :), Ω₀, Qj, λj, NSTEPS)
    end
    return ρmat
end

# compute an upper bound on ρ(eᵢ, Xₖ)
# ρmat is a matrix of size 2n x NSTEPS
# P is the eigenvectors matrix
# k is the time index

# `ρmat` is a matrix of size `2n × NSTEPS` such that
# the i-th row of `ρmat` contains all support function evaluations
# ρ(Φ^k dᵢ, Ω₀) for each k = 0, …, NSTEPS-1,
# where dᵢ = vᵢ and vᵢ is the i-th eigenvector of Φ^T if 1 ≤ i ≤ n and for i > n
# dᵢ = -vᵢ is the negative of the i-th eigenvector of Φ^T.
#
# This function receives:
# - `ρmat` -- matrix of support functions
# - `i`    -- an index `i` between `1` and `n`
# - `Qinv` -- the inverse of the matrix of eigenvectors of `Φ^T`, `Qinv`
# - `k`    -- time index `k` such that 1 ≤ k ≤ size(ρmat, 2) = NSTEPS
#
# The function returns an upper bound on ρ((Φ^T)^k eᵢ, Ω₀).
function _upper_bound_eig_dir(ρmat, i, Qinv, k)
    res = zero(eltype(Qinv))
    α = view(Qinv, :, i) # i-th column of Q^{-1}
    n = size(Qinv, 1)
    for (j, αⱼ) in enumerate(α)
        if αⱼ > 0
            res += αⱼ * ρmat[j, k]
        else
            res += abs(αⱼ) * ρmat[j + n, k]
        end
    end
    return res
end

# similar to _upper_bound_eig_dir but with eᵢ and -eᵢ
function _upper_bound_eig_dir_posneg(ρmat, i, Qinv, k)
    N = eltype(Qinv)
    a = zero(N)
    b = zero(N)

    α = view(Qinv, :, i) # i-th column of Q^{-1}
    n = size(Qinv, 1)
    for (j, αⱼ) in enumerate(α)
        if αⱼ > 0
            a += αⱼ * ρmat[j, k]
            b += αⱼ * ρmat[j + n, k]
        else
            a += abs(αⱼ) * ρmat[j + n, k]
            b += abs(αⱼ) * ρmat[j, k]
        end
    end
    return a, b
end

function _upper_bound_eig(ρmat::Matrix{N}, Qinv, NSTEPS) where {N}
    n = size(Qinv, 1)
    ρmat_box = Matrix{N}(undef, 2n, NSTEPS)
    @inbounds begin
        for k in 1:NSTEPS
            for i in 1:n
                a, b = _upper_bound_eig_dir_posneg(ρmat, i, Qinv, k)
                ρmat_box[i, k] = a
                ρmat_box[i + n, k] = b
            end
        end
    end
    return ρmat_box
end

# let x' = Ax, s.t. A has real eigenvalues of size n x n
# Aδ = P Λ P^{-1}, and Q = (P^{-1})^T
# Q has the eigenvectors of Φ^T = exp(A^T δ)
# the function assumes that the matrix Q is orthogonal and returns the box overapproximation
# of the flowipe by support function evaluation matrix
function reach_homog_eig_LGG09_box(Λ::Vector{N}, Q, Ω₀, NSTEPS) where {N}
    # both M₊ and M₋ are matrices of size n x NSTEPS
    # compute support functions along +vj where vj is column of Q (eigenvector of Φ^T)
    M₊ = reach_homog_eig_LGG09(Λ, Q, Ω₀, NSTEPS)

    # compute support functions along -vj where vj is column of Q (eigenvector of Φ^T)
    M₋ = reach_homog_eig_LGG09(Λ, -Q, Ω₀, NSTEPS)

    # center
    C = Q * (M₊ - M₋) / 2

    # radius
    R = abs.(Q) * (M₊ + M₋) / 2

    # support function of the hyperrectangular approximations
    B₊ = R + C
    B₋ = R - C

    return B₊, B₋
end

# ===================================
# Homogeneous case with invariant
# ===================================

function reach_homog_LGG09!(F::Vector{RT},
                            dirs::TN,
                            Ω₀::LazySet{N},
                            Φ::AbstractMatrix{N},
                            NSTEPS::Integer,
                            δ::N,
                            X::LazySet,
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

    # fill template reach-set sequence
    Δt = (zero(N) .. δ) + Δt0
    k = 1
    @inbounds while k <= NSTEPS

        # loop dispatch in threaded
        _reach_homog_LGG09_invariant!(Ω₀, rᵢ, ρmat, k, ndirs, threaded)

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

function _reach_homog_LGG09_invariant!(Ω₀, rᵢ, ρmat, k, ndirs, threaded::Val{false})
    for j in 1:ndirs
        d = view(rᵢ, :, j)
        ρmat[j, k] = ρ(d, Ω₀)
    end
end

function _reach_homog_LGG09_invariant!(Ω₀, rᵢ, ρmat, k, ndirs, threaded::Val{true})
    Threads.@threads for j in 1:ndirs
        d = view(rᵢ, :, j)
        ρmat[j, k] = ρ(d, Ω₀)
    end
end
