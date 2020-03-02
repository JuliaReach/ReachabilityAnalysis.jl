# ================
# Homogeneous case
# ================

using LinearAlgebra
# TODO: add to MathematicalSystems
Base.:*(im::IdentityMultiple, d::Diagonal) = im.M.λ * d

# X is the universal set => it is ignored
function reach_homog!(F::Vector{ReachSet{N, Zonotope{N}}},
                      Ω0::Zonotope{N},
                      Φ::AbstractMatrix,
                      NSTEPS::Integer,
                      δ::Float64,
                      max_order::Integer,
                      X::Universe) where {N}
    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    k = 2
    while k <= NSTEPS
        Rₖ = linear_map(Φ, set(F[k-1]))
        Rₖ = reduce_order(Rₖ, max_order)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end

# early termination checking for intersection with the invariant << TODO
function reach_homog!(F::Vector{ReachSet{N, Zonotope{N}}},
                      Ω0::Zonotope{N},
                      Φ::AbstractMatrix,
                      NSTEPS::Integer,
                      δ::Float64,
                      max_order::Integer,
                      X::LazySet) where {N}
    # initial reach set
    Δt = zero(N) .. δ
    F[1] = ReachSet(Ω0, Δt)

    k = 2
    while k <= NSTEPS
        Rₖ = linear_map(Φ, set(F[k-1]))
        Rₖ = reduce_order(Rₖ, max_order)
        Δt += δ
        F[k] = ReachSet(Rₖ, Δt)
        k += 1
    end
    return F
end

# ==================
# Inhomogeneous case
# ==================

# TODO: add invariant information

function reach_inhomog!(F::Vector{ReachSet{N, Zonotope{N}}},
                        Ω0::Zonotope{N},
                        Φ::AbstractMatrix,
                        NSTEPS::Integer,
                        δ::Float64,
                        max_order::Integer,
                        U::LazySet) where {N}
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
