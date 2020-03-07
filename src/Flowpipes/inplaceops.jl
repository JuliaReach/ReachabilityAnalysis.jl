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

# more efficient concrete projection for zonotopic sets (see LazySets#2013)
function LazySets.project(Z::AbstractZonotope{N}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    M = projection_matrix(collect(vars), dim(Z), N)
    return linear_map(M, Z)
    #return _project!(copy(Z), vars)
end

function project!(Z::AbstractZonotope{N}, vars::AbstractVector{Int}) where {N}
    M = projection_matrix(vars, dim(Z), N)
    return linear_map!(M, Z)
end

function project!(X::CartesianProduct{N, <:Interval, <:AbstractZonotope}, vars::AbstractVector{Int}) where {N}
    Z = convert(Zonotope, X)
    return project(Z, vars)
end

# fallback concrete projection
function _project(X::LazySet, vars::NTuple{D, T}) where {D, T<:Integer}
    return project(X, collect(vars))
end

# not available yet
function _project!(X::LazySet, vars::NTuple{D, T}) where {D, T<:Integer}
    error("the concrete inplace projection for a set of type $(typeof(X)) is not implemented yet")
#   return project!(X, collect(vars))
end

# more efficient concrete projection for zonotopic sets (see LazySets#2013)
function _project!(Z::AbstractZonotope{N}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    M = projection_matrix(collect(vars), dim(Z), N)
    return linear_map!(M, Z)
end

#=
# cf. https://github.com/JuliaArrays/ElasticArrays.jl
function minkowski_sum!(Zout::Zonotope{N, Vector{N}, ElasticArray{N}}, Z1::AbstractZonotope, Z2::AbstractZonotope) where {N}
    cnew = center(Zout)
    Gnew = genmat(Zout)
    cnew .= center(Z1) .+ center(Z2)
    append!(Gnew, genmat(Z1))
    append!(Gnew, genmat(Z2))
    Gnew .= hcat(genmat(Z1), genmat(Z2)) # TODO check if its: hcat(Gnew, ...?)
    return Zout
end
=#

#function reduce!(Z::AbstractZonotope, ord::Integer)
# ...
#end
