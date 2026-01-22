# ==================
# Inhomogeneous case
# ==================

# no invariant
function reach_inhomog_GLGM06!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                               Ω0::Zonotope{N,VN,MN},
                               Φ::AbstractMatrix,
                               NSTEPS::Integer,
                               δ::Float64,
                               max_order::Integer,
                               ::Universe,
                               U::LazySet,
                               reduction_method::AbstractReductionMethod,
                               Δt0::TimeInterval,
                               disjointness_method::AbstractDisjointnessMethod) where {N,VN,MN}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    Wk₊ = U
    Φ_power_k = similar(Φ)
    copyto!(Φ_power_k, Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = linear_map_minkowski_sum(Φ_power_k, Ω0, Wk₊)
        Rₖ = reduce_order(Rₖ, max_order, reduction_method)

        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        Wk₊ = linear_map_minkowski_sum(Φ_power_k, U, Wk₊)
        Wk₊ = reduce_order(Wk₊, max_order, reduction_method)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        Φ_power_k, Φ_power_k_cache = Φ_power_k_cache, Φ_power_k  # swap pointers
        k += 1
    end
    return F
end

# with invariant
function reach_inhomog_GLGM06!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                               Ω0::Zonotope{N,VN,MN},
                               Φ::AbstractMatrix,
                               NSTEPS::Integer,
                               δ::Float64,
                               max_order::Integer,
                               X::LazySet,
                               U::LazySet,
                               reduction_method::AbstractReductionMethod,
                               Δt0::TimeInterval,
                               disjointness_method::AbstractDisjointnessMethod) where {N,VN,MN}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    Wk₊ = U
    Φ_power_k = similar(Φ)
    copyto!(Φ_power_k, Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    @inbounds while k <= NSTEPS
        Rₖ = linear_map_minkowski_sum(Φ_power_k, Ω0, Wk₊)
        Rₖ = reduce_order(Rₖ, max_order, reduction_method)
        _isdisjoint(X, Rₖ, disjointness_method) && break
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        Wk₊ = linear_map_minkowski_sum(Φ_power_k, U, Wk₊)
        Wk₊ = reduce_order(Wk₊, max_order, reduction_method)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        Φ_power_k, Φ_power_k_cache = Φ_power_k_cache, Φ_power_k
        k += 1
    end
    if k < NSTEPS + 1
        resize!(F, k - 1)
    end
    return F
end
