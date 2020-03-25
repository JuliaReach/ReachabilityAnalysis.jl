# ================
# Homogeneous case
# ================

# X is the universal set => it is ignored
# in this case order reduction only has to be applied to the
# initial set because the loop does not create new generators
# it is assumed that order(Ω0) <= max_order
function reach_homog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                             Ω0::Zonotope{N, VN, MN},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::Float64,
                             max_order::Integer,
                             X::Universe) where {N, VN, MN}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = linear_map(Φ, set(F[k-1]))
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end

# homogeneous case using StaticArrays
# it is assumed that the order of Ω0 is at most max_order
function reach_homog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                             Ω0::Zonotope{N, VN, MN},
                             Φ::SMatrix,
                             NSTEPS::Integer,
                             δ::Float64,
                             max_order::Integer,
                             X::Universe) where {N, VN<:SVector, MN<:SMatrix}

    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    G0 = genmat(Ω0)
    c0 = center(Ω0)
    n, p = size(G0)
    np = n * p

    # cache for powers of Φ
    Φ_power_k = MMatrix{n, n, N, n * n}(Φ)
    Φ_power_k_cache = similar(Φ)

    k = 1
    @inbounds while k <= NSTEPS - 1
        ck = Φ_power_k * c0
        Gk = SMatrix{n, p, N, np}(Φ_power_k * G0)
        Zk = Zonotope(ck, Gk)
        Δt += δ
        k += 1
        F[k] = ReachSet(Zk, Δt)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
    end
    return F
end

# early termination checking for intersection with the invariant
# TODO : check stopping criterion
function reach_homog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                             Ω0::Zonotope{N, VN, MN},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::Float64,
                             max_order::Integer,
                             X::LazySet) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    k = 2
    while k <= NSTEPS
        Rₖ = linear_map(Φ, set(F[k-1]))
        is_intersection_empty(X, Rₖ) && break
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    if k < NSTEPS
        resize!(F, k)
    end
    return F
end

# TODO: static case with a state invariant


# OLD
#=
# in-place linear map
function reach_homog_GLGM06_inplace!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                                     Ω0::Zonotope{N, VN, MN},
                                     Φ::AbstractMatrix,
                                     NSTEPS::Integer,
                                     δ::Float64,
                                     max_order::Integer,
                                     X::Universe) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    Ω0red = reduce_order(Ω0, max_order)
    F[1] = ReachSet(Ω0red, Δt)
    c_Ω0 = center(Ω0red)
    G_Ω0 = genmat(Ω0red)
    n, p = size(G_Ω0)

    # cache for powers of Φ
    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    # vector of generator matrices
    cvec = Vector{VN}(undef, NSTEPS-1)
    Gvec = Vector{MN}(undef, NSTEPS-1)
    k = 1
    @inbounds while k <= NSTEPS-1
        cvec[k] = Vector{N}(undef, n)
        Gvec[k] = Matrix{N}(undef, n, p)
        mul!(cvec[k], Φ_power_k, c_Ω0)
        mul!(Gvec[k], Φ_power_k, G_Ω0)
        Rₖ = Zonotope(cvec[k], Gvec[k], remove_zero_generators=false)
        Δt += δ
        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
        F[k] = ReachSet(Rₖ, Δt)
    end
    return F
end
=#

#=
# in-place linear map
function reach_homog_GLGM06_inplace_v2!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                                     Ω0::Zonotope{N, VN, MN},
                                     Φ::AbstractMatrix,
                                     NSTEPS::Integer,
                                     δ::Float64,
                                     max_order::Integer,
                                     X::Universe) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    Ω0red = reduce_order(Ω0, max_order) # may not preserve SArrays

    c_Ω0 = center(Ω0red)
    G_Ω0 = genmat(Ω0red)
    n, p = size(G_Ω0)

    c_Ω0 = MVector{n}(c_Ω0)
    G_Ω0 = MMatrix{n, p}(G_Ω0)

    F[1] = ReachSet(Ω0red, Δt)

    k = 1
    @inbounds while k <= NSTEPS - 1
        #c_Rₖ = Φ * F[k].X.center
        mul!(F[k+1].X.center, Φ, F[k].X.center)
        mul!(F[k+1].X.generators, Φ, F[k].X.generators)
        #F[k+1].Δt
        #G_Rₖ = Φ * F[k].X.generators
        #Rₖ = Zonotope(c_Rₖ, G_Rₖ, remove_zero_generators=false)
        Δt += δ
        k += 1
        F[k] = ReachSet(Rₖ, Δt)
    end
    return F
end
=#

# ==================
# Inhomogeneous case
# ==================

# TODO: add case with invariant information
# TODO: check stopping criterion

function reach_inhomog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                               Ω0::Zonotope{N, VN, MN},
                               Φ::AbstractMatrix,
                               NSTEPS::Integer,
                               δ::Float64,
                               max_order::Integer,
                               X::LazySet,
                               U::LazySet) where {N, VN, MN}

    @warn "this function is still WIP"
    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    Wk₊ = U
    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    while k <= NSTEPS
        Rₖ = minkowski_sum(linear_map(Φ_power_k, Ω0), Wk₊, remove_zero_generators=false)
        Rₖ = reduce_order(Rₖ, max_order)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        Wk₊ = minkowski_sum(Wk₊, linear_map(Φ_power_k, U), remove_zero_generators=false)
        Wk₊ = reduce_order(Wk₊, max_order)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end
