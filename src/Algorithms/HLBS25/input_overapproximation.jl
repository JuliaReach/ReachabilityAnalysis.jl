function overapproximate_continuous_input(A::MatrixZonotope{N},
                                          B::MatrixZonotope{N},
                                          T::MatrixZonotope{N},
                                          U::SparsePolynomialZonotope{N},
                                          idg::IDGenerator,
                                          taylor_order::Int,
                                          A_norm::N;
                                          Δt::N) where {N}
    fresh!(idg, U)
    BU = overapproximate(B * U, SparsePolynomialZonotope)
    BUT = overapproximate(T * BU, SparsePolynomialZonotope)
    AT = A * T
    AΔt_norm = A_norm * Δt

    poly = LazySets.Approximations._truncated_operator_series(AT, BUT, taylor_order)
    Z = overapproximate(BU, Zonotope)
    rem = scale!(Δt, LazySets.Approximations._operator_series_remainder(Z, AΔt_norm, taylor_order))

    return remove_redundant_generators(minkowski_sum(poly, rem))
end

function overapproximate_discrete_input(A::MatrixZonotope{N},
                                        B::MatrixZonotope{N},
                                        U::SparsePolynomialZonotope{N},
                                        idg::IDGenerator,
                                        taylor_order::Int,
                                        A_norm::N,
                                        t::N) where {N}

    fresh!(idg, U)
    BUt = scale!(t, overapproximate(B * U, SparsePolynomialZonotope))
    At = scale(t, A)
    At_norm = A_norm * t
    #return LazySets.Approximations._overapproximate_emz_generic(At, BUt, taylor_order, At_norm)

    # use two separate buckets
    tayexp = LazySets.Approximations._truncated_operator_series(At, BUt, taylor_order)
    BUt_Z = overapproximate(BUt, Zonotope)
    lagrem = LazySets.Approximations._operator_series_remainder(BUt_Z, At_norm, taylor_order)
    return MinkowskiSum(tayexp, lagrem)
end
