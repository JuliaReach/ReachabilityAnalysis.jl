# overapproximate particular solution - discrete time
function overapproximate_discrete_input(A::MatrixZonotope{N},
                                        B::MatrixZonotope{N},
                                        U::SparsePolynomialZonotope{N},
                                        idg::IDGenerator,
                                        taylor_order::Int,
                                        A_norm::N,
                                        t::N) where {N}

    fresh!(idg, U)
    BUt = scale(overapproximate( B * U, SPZ), t)
    At = scale(A, t)
    At_norm = A_norm * t
    return LazySets.Approximations._overapproximate_emz_generic(At, BUt, taylor_order, At_norm)
end

function overapproximate_continuous_input(A::MatrixZonotope{N},
                                        B::MatrixZonotope{N},
                                        T::MatrixZonotope{N},
                                        U::SparsePolynomialZonotope{N},
                                        idg::IDGenerator,
                                        taylor_order::Int,
                                        A_norm::N;
                                        Δt::N) where {N}

    fresh!(idg, U)
    BU = overapproximate( B * U, SPZ)
    BUT = overapproximate( T * BU, SPZ)
    AT = A * T
    AΔt_norm = A_norm * Δt
    
    poly = LazySets.Approximations._truncated_operator_series(AT, BUT, taylor_order)
    Z = BU isa AbstractZonotope ? BU : overapproximate(BU, Zonotope)
    rem = scale(Δt, _operator_series_remainder(Z, AΔt_norm, taylor_order))

    return remove_redundant_generators(minkowski_sum(poly, rem))
end

function reach_inhomog_HLBS25!(F::Vector{ReachSet{N,S}},
                               Ω0::S,
                               Φ::MatrixZonotope{N,MN},
                               B::MatrixZonotope{N,MN},
                               U::S,
                               NSTEPS::Integer,
                               δ::N,
                               taylor_order::Integer,
                               Φ_norm::N,
                               recursive::Val{true},
                               max_order::Integer,
                               reduction_method::AbstractReductionMethod,
                               Δt0::IA.Interval,
                               idg::IDGenerator) where {N,S<:SparsePolynomialZonotope{N},
                                                        MN<:AbstractMatrix{N}}
    
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0

    t = sup(Δt)
    Pt = overapproximate_discrete_input(Φ, B, U, idg, taylor_order, Φ_norm, t)
    
    Hⱼ  = Ω0.X
    Pτⱼ = Ω0.Y
    Ω0 = concretize(Ω0)
    @inbounds F[1] = ReachSet(Ω0, Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))

    j = 1
    @inbounds while j < NSTEPS

        # update H(τⱼ)
        em = ExponentialMap(expΦδ, Hⱼ)
        Hⱼ₊₁ = overapproximate(em, SparsePolynomialZonotope, taylor_order)
        Hⱼ₊₁ = reduce_order(Hⱼ₊₁, max_order, reduction_method)

        PΔt = fresh(PΔt, A, B)

        # update P(τⱼ)
        em = ExponentialMap(expΦδ, Pτⱼ)
        tmp = overapproximate(em, SparsePolynomialZonotope, taylor_order)
        tmp = exact_sum(tmp, PΔt)
        Pτⱼ₊₁ = reduce_order(Pτ, max_order, reduction_method) 
    
        # update reach set
        Zⱼ₊₁ = exact_sum(Hⱼ₊₁, Pτⱼ₊₁)

        j += 1
        Δt += δ
        Pτⱼ = Pτⱼ₊₁
        Hⱼ = Hⱼ₊₁
        F[j] = ReachSet(Zⱼ₊₁, Δt)
    end
    return F
end