import .CorrectionHullMatrixZonotopeModule: overapproximate_discrete_input_split

@inline function _split_spz_components(P::SparsePolynomialZonotope{N}) where {N}
    n = dim(P)
    P_poly = SparsePolynomialZonotope(center(P), genmat_dep(P), zeros(N, n, 0), expmat(P),
                                      indexvector(P))
    P_zono = Zonotope(zeros(N, n), genmat_indep(P))
    return P_poly, P_zono
end

@inline function _merge_spz_components(P_poly::SparsePolynomialZonotope{N},
                                       P_zono::AbstractZonotope{N}) where {N}
    c = center(P_poly) + center(P_zono)
    P = SparsePolynomialZonotope(c, genmat_dep(P_poly), genmat(P_zono), expmat(P_poly),
                                 indexvector(P_poly))
    return remove_redundant_generators(P)
end

@inline function _reduce_buckets(P_poly::SparsePolynomialZonotope,
                                 P_zono::AbstractZonotope,
                                 max_order_poly::Integer,
                                 max_order_zono::Integer,
                                 reduction_method::AbstractReductionMethod)
    P_poly_r = reduce_order(P_poly, max_order_poly, reduction_method)
    P_zono_r = reduce_order(P_zono, max_order_zono, reduction_method)
    return P_poly_r, P_zono_r
end

@inline function _split_reduce_spz(P::SparsePolynomialZonotope,
                                   max_order_poly::Integer,
                                   max_order_zono::Integer,
                                   reduction_method::AbstractReductionMethod)
    P_poly, P_zono = _split_spz_components(P)
    return _reduce_buckets(P_poly, P_zono, max_order_poly, max_order_zono, reduction_method)
end

function _prepare_discrete_input_buckets(Φ::MatrixZonotope{N,MN},
                                         B::MatrixZonotope{N,MN},
                                         U::SparsePolynomialZonotope{N},
                                         idg::IDGenerator,
                                         taylor_order::Integer,
                                         Φ_norm::N,
                                         t::N,
                                         max_order_poly::Integer,
                                         max_order_zono::Integer,
                                         reduction_method::AbstractReductionMethod) where {N,
                                                                                           MN<:AbstractMatrix{N}}
    PΔt_tay, PΔt_rem = overapproximate_discrete_input_split(Φ, B, U, idg, taylor_order, Φ_norm, t)
    PΔt_poly, PΔt_zono_tay = _split_spz_components(PΔt_tay)
    PΔt_zono = remove_redundant_generators(minkowski_sum(PΔt_zono_tay, PΔt_rem))
    return _reduce_buckets(PΔt_poly, PΔt_zono, max_order_poly, max_order_zono, reduction_method)
end
