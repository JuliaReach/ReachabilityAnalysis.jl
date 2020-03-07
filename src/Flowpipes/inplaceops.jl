# Extension of some common LazySets operations, but without new allocations.

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

function project!(Z::AbstractZonotope{N}, vars::AbstractVector{Int}) where {N}
    M = projection_matrix(vars, dim(Z), N)
    return linear_map!(M, Z)
end

function project!(X::CartesianProduct{N, <:Interval, <:AbstractZonotope}, vars::AbstractVector{Int}) where {N}
    Z = convert(Zonotope, X)
    return project(Z, vars)
end

function minkowski_sum!(Zout::AbstractZonotope, Z1::AbstractZonotope, Z2::AbstractZonotope)
    cnew = center(Zout)
    Gnew = genmat(Zout)
    cnew .= center(Z1) .+ center(Z2)
    Gnew .= hcat(genmat(Z1), genmat(Z2)) # TODO check if its: hcat(Gnew, ...?)
    return Zout
end

#function reduce!(Z::AbstractZonotope, ord::Integer)
# ...
#end
