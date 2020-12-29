# ===========================
# Homogeneous case
# ===========================

function reach_homog_VREP!(F::Vector{ReachSet{N, VP}},
                           Ω0::VP,
                           Φ::AbstractMatrix,
                           NSTEPS::Integer,
                           δ::Float64,
                           X::Universe,
                           Δt0::TimeInterval) where {N, VN, VP<:VPOLY{N, VN}}

    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
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
