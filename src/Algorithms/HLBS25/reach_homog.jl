function reach_homog_HLBS25!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                             Ω0::Zonotope{N,VN,MN},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::N,
                             max_order::Integer,
                             taylor_order::Integer,
                             recursive::Val{false},
                             reduction_method::AbstractReductionMethod,
                             Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))
    expΦδ = overapproximate(expΦδ, MatrixZonotope, taylor_order)

    j = 2
    @inbounds while j ≤ NSTEPS
        # current set
        Zⱼ = set(F[j - 1])
        Zⱼ₊₁ = overapproximate(expΦδ * Zⱼ, Zonotope)
        Zⱼ₊₁ʳ = reduce_order(Zⱼ₊₁, max_order, reduction_method)

        j += 1
        Δt += δ
        F[j] = ReachSet(Zⱼ₊₁ʳ, Δt)
    end
    return F
end

function reach_homog_HLBS25!(F::Vector{ReachSet{N,Zonotope{N,VN,MN}}},
                             Ω0::Zonotope{N,VN,MN},
                             Φ::AbstractMatrix,
                             NSTEPS::Integer,
                             δ::N,
                             max_order::Integer,
                             taylor_order::Integer,
                             recursive::Val{true},
                             reduction_method::AbstractReductionMethod,
                             Δt0::TN) where {N,TN,VN,MN}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    j = 2
    @inbounds while j ≤ NSTEPS
        Zⱼ = set(F[j - 1])

        expΦδ = MatrixZonotopeExp(scale(δ, Φ))
        em = ExponentialMap(expΦδ, Zⱼ)

        Zⱼ₊₁ = overapproximate(em, Zonotope, taylor_order)
        Zⱼ₊₁ʳ = reduce_order(Zⱼ₊₁, max_order, reduction_method)

        j += 1
        Δt += δ
        F[j] = ReachSet(Zⱼ₊₁ʳ, Δt)
    end
    return F
end
