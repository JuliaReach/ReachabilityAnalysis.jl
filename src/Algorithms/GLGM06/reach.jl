# ================
# Homogeneous case
# ================

function reach_homog!(R::Vector{ReachSet{Zonotope{Float64}}},
                      Ω0::Zonotope,
                      Φ::AbstractMatrix,
                      N::Int,
                      δ::Float64,
                      max_order::Int)
    # initial reach set
    0,
    δ
    R[1] = ReachSet(Ω0, t0..t1)

    k = 2
    while k <= N
        Rₖ = linear_map(Φ, set(R[k-1]))
        Rₖ = reduce_order(Rₖ, max_order)
        t0 = t1
        t1 += δ
        R[k] = ReachSet(Rₖ, t0..t1)

        k += 1
    end
    return R
end

# ==================
# Inhomogeneous case
# ==================

function reach_inhomog!(R::Vector{ReachSet{Zonotope{Float64}}},
                        Ω0::Zonotope,
                        U::LazySet,
                        Φ::AbstractMatrix,
                        N::Int,
                        δ::Float64,
                        max_order::Int)
    # initial reach set
    t0, t1 = zero(δ), δ
    R[1] = ReachSet(Ω0, t0, t1)

    Wk₊ = U
    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    k = 2
    while k <= N
        Rₖ = minkowski_sum(linear_map(Φ_power_k, Ω0), Wk₊)
        Rₖ = reduce_order(Rₖ, max_order)
        t0 = t1
        t1 += δ
        R[k] = ReachSet(Rₖ, t0, t1)

        Wk₊ = minkowski_sum(Wk₊, linear_map(Φ_power_k, U))
        Wk₊ = reduce_order(Wk₊, max_order)

        mul!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end
    return R
end
