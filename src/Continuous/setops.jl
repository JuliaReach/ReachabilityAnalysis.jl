# compute `(M * Z1) ⊕ Z2`
function linear_map_minkowski_sum(M::AbstractMatrix, X::LazySet, Y::LazySet)
    return minkowski_sum(linear_map(M, X), Y)
end

# no intermediate allocations
function linear_map_minkowski_sum(M::Matrix, Z1::AbstractZonotope, Z2::AbstractZonotope)
    m1, n1 = size(M)
    G1 = genmat(Z1)
    m2, n2 = size(G1)
    G2 = genmat(Z2)
    m3, n3 = size(G2)
    @assert n1 == m2 && m1 == m3 "invalid dimensions: $(size(M)), $(size(Z1)), $(size(Z2))"
    N = promote_type(eltype(M), eltype(Z1), eltype(Z2))

    c = Vector{N}(undef, m1)
    mul!(c, M, center(Z1))
    c .+= center(Z2)

    G = Matrix{N}(undef, m1, n2 + n3)
    mul!(view(G, :, 1:n2), M, G1)
    G[:, (n2 + 1):end] = G2

    return Zonotope(c, G)
end
