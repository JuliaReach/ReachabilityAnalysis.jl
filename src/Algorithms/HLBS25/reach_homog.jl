function reach_homog_HLBS25!(F::Vector{ReachSet{N,S}},
                             Ω0::S,
                             Φ::MatrixZonotope{N,MN},
                             NSTEPS::Integer,
                             δ::N,
                             taylor_order::Integer,
                             recursive::Val{false},
                             max_order_poly::Integer,
                             max_order_zono::Integer,
                             reduction_method::AbstractReductionMethod,
                             Δt0::TimeInterval,
                             idg::IDGenerator) where {N,S<:SparsePolynomialZonotope{N},
                                                      MN<:AbstractMatrix{N}}
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Z_poly, Z_zono = _split_reduce_spz(Ω0, max_order_poly, max_order_zono, reduction_method)
    @inbounds F[1] = ReachSet(_merge_spz_components(Z_poly, Z_zono), Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))
    expΦδ_approx = overapproximate(expΦδ, MatrixZonotope, taylor_order)
    synchronize!(idg, Φ, Z_poly)
    fresh!(idg, expΦδ_approx)

    j = 1
    @inbounds while j < NSTEPS
        Z_poly_next = overapproximate(expΦδ_approx * Z_poly, SparsePolynomialZonotope)
        Z_zono_next = overapproximate(expΦδ_approx * Z_zono, Zonotope)
        Z_poly_next, Z_zono_next = _reduce_buckets(Z_poly_next, Z_zono_next, max_order_poly,
                                                   max_order_zono, reduction_method)

        j += 1
        Δt += δ
        Z_poly = Z_poly_next
        Z_zono = Z_zono_next
        F[j] = ReachSet(_merge_spz_components(Z_poly, Z_zono), Δt)
    end
    return F
end

function reach_homog_HLBS25!(F::Vector{ReachSet{N,S}},
                             Ω0::S,
                             Φ::MatrixZonotope{N,MN},
                             NSTEPS::Integer,
                             δ::N,
                             taylor_order::Integer,
                             recursive::Val{true},
                             max_order_poly::Integer,
                             max_order_zono::Integer,
                             reduction_method::AbstractReductionMethod,
                             Δt0::TimeInterval,
                             idg::IDGenerator) where {N,S<:SparsePolynomialZonotope{N},
                                                      MN<:AbstractMatrix{N}}
    _ = idg
    # initial reach set
    Δt = (zero(N) .. δ) + Δt0
    Z_poly, Z_zono = _split_reduce_spz(Ω0, max_order_poly, max_order_zono, reduction_method)
    @inbounds F[1] = ReachSet(_merge_spz_components(Z_poly, Z_zono), Δt)

    expΦδ = MatrixZonotopeExp(scale(δ, Φ))

    j = 1
    @inbounds while j < NSTEPS
        Z_poly_next = overapproximate(ExponentialMap(expΦδ, Z_poly), SparsePolynomialZonotope,
                                      taylor_order)
        Z_zono_next = overapproximate(ExponentialMap(expΦδ, Z_zono), Zonotope, taylor_order)
        Z_poly_next, Z_zono_next = _reduce_buckets(Z_poly_next, Z_zono_next, max_order_poly,
                                                   max_order_zono, reduction_method)

        j += 1
        Δt += δ
        Z_poly = Z_poly_next
        Z_zono = Z_zono_next
        F[j] = ReachSet(_merge_spz_components(Z_poly, Z_zono), Δt)
    end
    return F
end
