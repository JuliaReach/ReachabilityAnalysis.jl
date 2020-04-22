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
                               U::LazySet,
                               reduction_method::AbstractReductionMethod) where {N, VN, MN}

    @warn "this function is still WIP"
    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    Wk₊ = U
    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    while k <= NSTEPS
        Rₖ = _minkowski_sum(_linear_map(Φ_power_k, Ω0), Wk₊)
        Rₖ = _reduce_order(Rₖ, max_order, reduction_method)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)

        Wk₊ = _minkowski_sum(Wk₊, _linear_map(Φ_power_k, U))
        Wk₊ = _reduce_order(Wk₊, max_order, reduction_method)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return F
end
