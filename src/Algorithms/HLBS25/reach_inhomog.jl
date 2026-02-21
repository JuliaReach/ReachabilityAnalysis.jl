function reach_inhomog_HLBS25!(F::Vector{ReachSet{N,ZS}},
                               Ω0::ExactSum{N,S1,S2},
                               Φ::MatrixZonotope{N,MN},
                               B::MatrixZonotope{N,MN},
                               U::SparsePolynomialZonotope{N},
                               NSTEPS::Integer,
                               δ::N,
                               taylor_order::Integer,
                               Φ_norm::N,
                               recursive::Val{true},
                               max_order::Integer,
                               reduction_method::AbstractReductionMethod,
                               Δt0::IA.Interval,
                               idg::IDGenerator) where {N,ZS<:SparsePolynomialZonotope{N},
                                                        S1<:SparsePolynomialZonotope{N},
                                                        S2<:SparsePolynomialZonotope{N},
                                                        MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(concretize(Ω0), Δt)

    t = sup(Δt)
    PΔt = overapproximate_discrete_input(Φ, B, U, idg, taylor_order, Φ_norm, t)
    PΔt = concretize(PΔt)

    Hⱼ = Ω0.X
    Pτⱼ = Ω0.Y

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))

    j = 1
    @inbounds while j < NSTEPS
        # update H(τⱼ)
        emH = ExponentialMap(expΦδ, Hⱼ)
        Hⱼ₊₁ = overapproximate(emH, SparsePolynomialZonotope, taylor_order)
        Hⱼ₊₁ = reduce_order(Hⱼ₊₁, max_order, reduction_method)

        # update P(τⱼ)
        emP = ExponentialMap(expΦδ, Pτⱼ)
        Ptmp = overapproximate(emP, SparsePolynomialZonotope, taylor_order)
        Pτⱼ₊₁ = exact_sum(Ptmp, PΔt)
        Pτⱼ₊₁ = overapproximate(Pτⱼ₊₁, SparsePolynomialZonotope)
        Pτⱼ₊₁ = reduce_order(Pτⱼ₊₁, max_order, reduction_method)

        # update reach set as concrete SPZ (ExactSum is only an internal holder)
        Zⱼ₊₁ = exact_sum(Hⱼ₊₁, Pτⱼ₊₁)

        j += 1
        Δt += δ
        Pτⱼ = Pτⱼ₊₁
        Hⱼ = Hⱼ₊₁
        F[j] = ReachSet(Zⱼ₊₁, Δt)
    end
    return F
end

function reach_inhomog_HLBS25!(F::Vector{ReachSet{N,ZS}},
                               Ω0::ExactSum{N,S1,S2},
                               Φ::MatrixZonotope{N,MN},
                               B::MatrixZonotope{N,MN},
                               U::SparsePolynomialZonotope{N},
                               NSTEPS::Integer,
                               δ::N,
                               taylor_order::Integer,
                               Φ_norm::N,
                               recursive::Val{false},
                               max_order::Integer,
                               reduction_method::AbstractReductionMethod,
                               Δt0::IA.Interval,
                               idg::IDGenerator) where {N,ZS<:SparsePolynomialZonotope{N},
                                                        S1<:SparsePolynomialZonotope{N},
                                                        S2<:SparsePolynomialZonotope{N},
                                                        MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    @inbounds F[1] = ReachSet(concretize(Ω0), Δt)

    t = sup(Δt)
    PΔt = overapproximate_discrete_input(Φ, B, U, idg, taylor_order, Φ_norm, t)

    Hⱼ = Ω0.X
    Pτⱼ = Ω0.Y

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))
    expΦδ_approx = overapproximate(expΦδ, MatrixZonotope, taylor_order)
    j = 1
    @inbounds while j < NSTEPS
        # update H(τⱼ)
        Hⱼ₊₁ = overapproximate(expΦδ_approx * Hⱼ, SparsePolynomialZonotope)
        Hⱼ₊₁ = reduce_order(Hⱼ₊₁, max_order, reduction_method)

        # update P(τⱼ)
        Ptmp = overapproximate(expΦδ_approx * Pτⱼ, SparsePolynomialZonotope)
        Pτⱼ₊₁ = exact_sum(Ptmp, PΔt)
        Pτⱼ₊₁ = overapproximate(Pτⱼ₊₁, SparsePolynomialZonotope)
        Pτⱼ₊₁ = reduce_order(Pτⱼ₊₁, max_order, reduction_method)

        # update reach set as concrete SPZ (ExactSum is only an internal holder)
        Zⱼ₊₁ = exact_sum(Hⱼ₊₁, Pτⱼ₊₁)

        j += 1
        Δt += δ
        Pτⱼ = Pτⱼ₊₁
        Hⱼ = Hⱼ₊₁
        F[j] = ReachSet(Zⱼ₊₁, Δt)
    end
    return F
end
