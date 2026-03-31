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
                               max_order_poly::Integer,
                               max_order_zono::Integer,
                               reduction_method::AbstractReductionMethod,
                               Δt0::IA.Interval,
                               idg::IDGenerator) where {N,ZS<:SparsePolynomialZonotope{N},
                                                        S1<:SparsePolynomialZonotope{N},
                                                        S2<:SparsePolynomialZonotope{N},
                                                        MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    t = sup(Δt)

    H_poly, H_zono = _split_spz_components(Ω0.X)
    Pτ_poly, Pτ_zono = _split_spz_components(Ω0.Y)

    @inbounds F[1] = ReachSet(concretize(Ω0), Δt)

    PΔt_poly, PΔt_zono = _prepare_discrete_input_buckets(Φ, B, U, idg, taylor_order, Φ_norm, t,
                                                         max_order_poly, max_order_zono,
                                                         reduction_method)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))

    j = 1
    @inbounds while j < NSTEPS

        if j%10 == 0
            println("Step $j / $NSTEPS, time = $(round(sup(Δt), digits=3))")
        end
        # update H(τⱼ) with separate polynomial and zonotopic buckets
        H_poly_next = overapproximate(ExponentialMap(expΦδ, H_poly), SparsePolynomialZonotope,
                                      taylor_order)
        H_zono_next = overapproximate(ExponentialMap(expΦδ, H_zono), Zonotope, taylor_order)
        H_poly_next, H_zono_next = _reduce_buckets(H_poly_next, H_zono_next, max_order_poly,
                                                   max_order_zono, reduction_method)

        # update P(τⱼ) with separate polynomial and zonotopic buckets
        # Preserve the algorithm's FRESH(P(Δt), A, B) step on the polynomial bucket.
        fresh!(idg, PΔt_poly, Φ, B)
        Pτ_poly_tmp = overapproximate(ExponentialMap(expΦδ, Pτ_poly), SparsePolynomialZonotope,
                                      taylor_order)
        Pτ_zono_tmp = overapproximate(ExponentialMap(expΦδ, Pτ_zono), Zonotope, taylor_order)
        Pτ_poly_next = exact_sum(Pτ_poly_tmp, PΔt_poly)
        Pτ_zono_next = minkowski_sum(Pτ_zono_tmp, PΔt_zono)
        Pτ_poly_next, Pτ_zono_next = _reduce_buckets(Pτ_poly_next, Pτ_zono_next, max_order_poly,
                                                     max_order_zono, reduction_method)

        # merge only for saving the reach-set
        H_set = _merge_spz_components(H_poly_next, H_zono_next)
        Pτ_set = _merge_spz_components(Pτ_poly_next, Pτ_zono_next)
        Zⱼ₊₁ = remove_redundant_generators(exact_sum(H_set, Pτ_set))

        j += 1
        Δt += δ
        H_poly = H_poly_next
        H_zono = H_zono_next
        Pτ_poly = Pτ_poly_next
        Pτ_zono = Pτ_zono_next
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
                               max_order_poly::Integer,
                               max_order_zono::Integer,
                               reduction_method::AbstractReductionMethod,
                               Δt0::IA.Interval,
                               idg::IDGenerator) where {N,ZS<:SparsePolynomialZonotope{N},
                                                        S1<:SparsePolynomialZonotope{N},
                                                        S2<:SparsePolynomialZonotope{N},
                                                        MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    t = sup(Δt)

    H_poly, H_zono = _split_spz_components(Ω0.X)
    Pτ_poly, Pτ_zono = _split_spz_components(Ω0.Y)

    @inbounds F[1] = ReachSet(concretize(Ω0), Δt)

    PΔt_poly, PΔt_zono = _prepare_discrete_input_buckets(Φ, B, U, idg, taylor_order, Φ_norm, t,
                                                         max_order_poly, max_order_zono,
                                                         reduction_method)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))
    expΦδ_approx = overapproximate(expΦδ, MatrixZonotope, taylor_order)
    synchronize!(idg, Φ, B, U, H_poly, Pτ_poly, PΔt_poly)
    fresh!(idg, expΦδ_approx)

    j = 1
    @inbounds while j < NSTEPS

        if j%10 == 0
            println("Step $j / $NSTEPS, time = $(round(sup(Δt), digits=3))")
        end
        # update H(τⱼ) with separate polynomial and zonotopic buckets
        H_poly_next = overapproximate(expΦδ_approx * H_poly, SparsePolynomialZonotope)
        H_zono_next = overapproximate(expΦδ_approx * H_zono, Zonotope)
        H_poly_next, H_zono_next = _reduce_buckets(H_poly_next, H_zono_next, max_order_poly,
                                                   max_order_zono, reduction_method)

        # update P(τⱼ) with separate polynomial and zonotopic buckets
        fresh!(idg, PΔt_poly, Φ, B)
        Pτ_poly_tmp = overapproximate(expΦδ_approx * Pτ_poly, SparsePolynomialZonotope)
        Pτ_zono_tmp = overapproximate(expΦδ_approx * Pτ_zono, Zonotope)
        Pτ_poly_next = exact_sum(Pτ_poly_tmp, PΔt_poly)
        Pτ_zono_next = minkowski_sum(Pτ_zono_tmp, PΔt_zono)
        Pτ_poly_next, Pτ_zono_next = _reduce_buckets(Pτ_poly_next, Pτ_zono_next, max_order_poly,
                                                     max_order_zono, reduction_method)

        # merge only for saving the reach-set
        H_set = _merge_spz_components(H_poly_next, H_zono_next)
        Pτ_set = _merge_spz_components(Pτ_poly_next, Pτ_zono_next)
        Zⱼ₊₁ = remove_redundant_generators(exact_sum(H_set, Pτ_set))

        j += 1
        Δt += δ
        H_poly = H_poly_next
        H_zono = H_zono_next
        Pτ_poly = Pτ_poly_next
        Pτ_zono = Pτ_zono_next
        F[j] = ReachSet(Zⱼ₊₁, Δt)
    end
    return F
end
