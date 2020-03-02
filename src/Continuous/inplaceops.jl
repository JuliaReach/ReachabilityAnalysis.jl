# common set operations, in-place

function scale!(α::Real, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)

    c .= α .* c
    G .= α .* G
    return Zonotope(c, G)
end

function linear_map!(M::AbstractMatrix, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)

    c .= M * c
    G .= M .* G
    return Z
end

#function reduce!(Z::AbstractZonotope, ord::Integer)
# ...
#end
