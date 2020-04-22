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
                             X::Universe,
                             preallocate::Val{false}) where {N, VN, MN}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = _linear_map(Φ, set(F[k-1]))
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end

# version that preallocates the output zonotopes
function reach_homog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}},
                             Ω0::Zonotope{N, Vector{N}, Matrix{N}},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::Float64,
                             max_order::Integer,
                             X::Universe,
                             preallocate::Val{true}) where {N}

    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    n, p = size(Ω0.generators)

    # preallocate output
    Zout = Vector{Zonotope{N, Vector{N}, Matrix{N}}}(undef, NSTEPS)

    if p == 0
        Zout[1] = Zonotope(Ω0.center, zeros(N, n, 1))
        p = 1
    else
        Zout[1] = Ω0
    end

    @inbounds for i in 2:NSTEPS
        c = Vector{N}(undef, n)
        G = Matrix{N}(undef, n, p)
        Zout[i] = Zonotope(c, G)
    end

    k = 2
    @inbounds while k <= NSTEPS
        _linear_map!(Zout[k], Φ, Zout[k-1])
        Δt += δ
        F[k] = ReachSet(Zout[k], Δt)
        k += 1
    end
    return F
end

# check interection with invariant on the loop
# TODO: add variation with `preallocate` option, so this can be used
function reach_homog_GLGM06!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                             Ω0::Zonotope{N, VN, MN},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::Float64,
                             max_order::Integer,
                             X::LazySet) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    @inbounds F[1] = ReachSet(Ω0, Δt)

    k = 2
    while k <= NSTEPS
        Rₖ = linear_map(Φ, set(F[k-1]))
        is_intersection_empty(X, Rₖ) && break
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k)
    end
    return F
end

#= O L D
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
=#

#=
# follows the implementation for the general case
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

    k = 2
    @inbounds while k <= NSTEPS
        # TODO: benchmark wrt other method
        # .. substitute with doing Rₖ = linear_map(Φ, set(F[k-1]))
        X = set(F[k-1])
        ck = Φ * X.center
        Gk = Φ * X.generators

        Δt += δ
        Zk = Zonotope(ck, Gk)
        F[k] = ReachSet(Zk, Δt)
        k += 1
    end
    return F
end
=#

# OLD
#=
# in-place linear map
function reach_homog_GLGM06_inplace!(F::Vector{ReachSet{N, Zonotope{N, VN, MN}}},
                                     Ω0::Zonotope{N, VN, MN},
                                     Φ::AbstractMatrix,
                                     NSTEPS::Integer,
                                     δ::Float64,
                                     max_order::Integer,
                                     X::Universe,
                                     reduction_method::AbstractReductionMethod) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    Ω0red = reduce_order(Ω0, max_order, reduction_method)
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
        Rₖ = Zonotope(cvec[k], Gvec[k])
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
                                     X::Universe,
                                     reduction_method::AbstractReductionMethod) where {N, VN, MN}
    # initial reach set
    Δt = zero(N) .. δ
    Ω0red = reduce_order(Ω0, max_order, reduction_method) # may not preserve SArrays

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
        #Rₖ = Zonotope(c_Rₖ, G_Rₖ)
        Δt += δ
        k += 1
        F[k] = ReachSet(Rₖ, Δt)
    end
    return F
end
=#
