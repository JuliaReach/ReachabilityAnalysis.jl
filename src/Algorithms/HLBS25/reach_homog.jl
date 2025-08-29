function reach_homog_HLBS25!(F::Vector{ReachSet{N,S}},
                             Ω0::S,
                             Φ::MatrixZonotope{N,MN},
                             NSTEPS::Integer,
                             δ::N,
                             max_order::Integer,
                             taylor_order::Integer,
                             ::Val{false},
                             reduction_method::AbstractReductionMethod,
                             Δt0::IA.Interval) where {N,S<:SparsePolynomialZonotope{N},
                                                      MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))
    expΦδ_approx = overapproximate(expΦδ, MatrixZonotope, taylor_order)

    j = 1
    @inbounds while j < NSTEPS
        Zⱼ = set(F[j])
        Zⱼ₊₁ = overapproximate(expΦδ_approx * Zⱼ, SparsePolynomialZonotope)
        Zⱼ₊₁ʳ = reduce_order(Zⱼ₊₁, max_order, reduction_method)

        j += 1
        Δt += δ
        F[j] = ReachSet(Zⱼ₊₁ʳ, Δt)
    end
    return F
end

function reach_homog_HLBS25!(F::Vector{ReachSet{N,S}},
                             Ω0::S,
                             Φ::MatrixZonotope{N,MN},
                             NSTEPS::Integer,
                             δ::N,
                             max_order::Integer,
                             taylor_order::Integer,
                             ::Val{true},
                             reduction_method::AbstractReductionMethod,
                             Δt0::IA.Interval) where {N,S<:SparsePolynomialZonotope{N},
                                                      MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(Ω0, Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))

    j = 1
    @inbounds while j < NSTEPS
        Zⱼ = set(F[j])

        em = ExponentialMap(expΦδ, Zⱼ)

        Zⱼ₊₁ = overapproximate(em, SparsePolynomialZonotope, taylor_order)
        Zⱼ₊₁ʳ = reduce_order(Zⱼ₊₁, max_order, reduction_method)

        j += 1
        Δt += δ
        F[j] = ReachSet(Zⱼ₊₁ʳ, Δt)
    end
    return F
end
