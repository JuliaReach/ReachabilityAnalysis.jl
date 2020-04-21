# =========================
# (Re) conversion
# =========================

# no-op
_reconvert(Ω0::Zonotope{N, Vector{N}, Matrix{N}}, static::Val{false}) where {N} = Ω0
_reconvert(Ω0::Zonotope{N, <:SVector, <:SMatrix}, static::Val{true}) where {N} = Ω0

# convert any zonotope to be represented wih regular arrays
function _reconvert(Ω0::Zonotope, static::Val{false})
    Ω0 = Zonotope(Vector(Ω0.center), Matrix(Ω0.generators))
end

# convert any zonotope to be represented with static arrays
function _reconvert(Ω0::Zonotope{N, VN, MN}, static::Val{true}) where {N, VN, MN}
    n, p = size(Ω0.generators) # dimension and number of generators
    Ω0 = Zonotope(SVector{n, N}(Ω0.center), SMatrix{n, p, N, n*p}(Ω0.generators))
end

# no-op
_reconvert(Ω0::Hyperrectangle{N, Vector{N}, Vector{N}}, static::Val{false}) where {N} = Ω0
_reconvert(Ω0::Hyperrectangle{N, <:SVector, <:SVector}, static::Val{true}) where {N} = Ω0

# convert any Hyperrectangle to be represented wih regular arrays
function _reconvert(Ω0::Hyperrectangle, static::Val{false})
    Ω0 = Hyperrectangle(Vector(Ω0.center), Matrix(Ω0.radius), check_bounds=false)
end

# convert any Hyperrectangle to be represented with static arrays
function _reconvert(Ω0::Hyperrectangle{N, VNC, VNR}, static::Val{true}) where {N, VNC, VNR}
    n = length(Ω0.center) # dimension
    Ω0 = Hyperrectangle(SVector{n, N}(Ω0.center), SVector{n, N}(Ω0.radius), check_bounds=false)
end

# no-op
_reconvert(Φ::Matrix{N}, static::Val{false}) where {N} = Φ
_reconvert(Φ::IntervalMatrix{N}, static::Val{false}) where {N} = Φ
_reconvert(Φ::AbstractMatrix, static::Val{false}) = Matrix(Φ)
_reconvert(Φ::SMatrix, static::Val{true}) = Φ

function _reconvert(Φ::AbstractMatrix{N}, static::Val{true}) where {N}
    n = size(Φ, 1)
    Φ = SMatrix{n, n, N, n*n}(Φ)
end

function _reconvert(Φ::IntervalMatrix{N, IN, Matrix{IN}}, static::Val{true}) where {N, IN}
    n = size(Φ, 1)
    Φ = IntervalMatrix(SMatrix{n, n, IN, n*n}(Φ))
end

# fallback implementation for conversion (if applicable) or overapproximation
function _convert_or_overapproximate(T::Type{<:AbstractPolytope}, X::LazySet)
    if applicable(convert, T, X)
        return convert(T, X)
    elseif applicable(overapproximate, X, T)
        return overapproximate(X, T)
    else
        return convert(T, overapproximate(X, Hyperrectangle))
    end
end

#=
# no-op
function _convert_or_overapproximate2(T::Type{<:Zonotope}, X::Zonotope)
     return X
end
=#
#=
function _convert_or_overapproximate2(T::Type{Zonotope}, X::Interval)
    return convert(T, X)
end
=#

function _convert_or_overapproximate(X::LazySet, T::Type{<:AbstractPolytope})
    return _convert_or_overapproximate(T, X)
end

# =========================
# In-place ops
# =========================

# Extension of some common LazySets operations, some of them in-place

function _minkowski_sum(Z1::Zonotope{N}, Z2::Zonotope{N}) where {N}
    cnew = center(Z1) + center(Z2)
    Gnew = hcat(genmat(Z1), genmat(Z2))
    return Zonotope(cnew, Gnew)
end

# in-place scale of a zonotope
function scale!(α::Real, Z::Zonotope)
    c = Z.center
    G = Z.generators
    c .= α .* c
    G .= α .* G
    return Z
end

# in-place linear map of a zonotope
@inline function _linear_map!(Zout::Zonotope, M::AbstractMatrix, Z::Zonotope)
    mul!(Zout.center, M, Z.center)
    mul!(Zout.generators, M, Z.generators)
    return Zout
end

@inline function _linear_map(M::AbstractMatrix, Z::Zonotope)
    cout = M * Z.center
    Gout = M * Z.generators
    return Zonotope(cout, Gout)
end

# =========================
# Projection
# =========================

# fallback concrete projection
function _project(X::LazySet, vars::NTuple{D, T}) where {D, T<:Integer}
    return LazySets.project(X, collect(vars))
end

# more efficient concrete projection for zonotopic sets (see LazySets#2013)
function _project(Z::AbstractZonotope{N}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    M = projection_matrix(collect(vars), dim(Z), N)
    return linear_map(M, Z)  # return _project!(copy(Z), vars)
end

# cartesian product of an interval and a hyperrectangular set
function _project(cp::CartesianProduct{N, Interval{N, IA.Interval{N}}, <:AbstractHyperrectangle{N}}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    # NOTE: can be optimized further, eg. for the case that 1 ∉ vars and
    # without copying the full center and radius vectors of H
    Δt = cp.X
    H = cp.Y
    cH = vcat(center(Δt), center(H))
    rH = vcat(radius_hyperrectangle(Δt), radius_hyperrectangle(H))
    vars_vec = collect(vars)
    return Hyperrectangle(cH[vars_vec], rH[vars_vec], check_bounds=false)
end

# hyperrectangular set
function _project(H::AbstractHyperrectangle{N}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    cH = center(H)
    rH = radius_hyperrectangle(H)
    vars_vec = collect(vars)
    return Hyperrectangle(cH[vars_vec], rH[vars_vec], check_bounds=false)
end

# cartesian product of an interval and a zonotopic set
function _project(cp::CartesianProduct{N, Interval{N, IA.Interval{N}}, <:AbstractZonotope{N}}, vars::NTuple{D, T}) where {N, D, T<:Integer}
    vars_vec = collect(vars)
    M = projection_matrix(vars_vec, dim(cp), N)
    Z = convert(Zonotope, cp)
    return linear_map(M, Z)
end

#=
# TODO: add to LazySets
# specialize to cartesian product of Interval vs Zonotopic set
function Base.convert(::Type{Zonotope}, cp::CartesianProduct{N, <:AbstractZonotope{N}, <:AbstractZonotope{N}}) where {N<:Real}
    Z1, Z2 = cp.X, cp.Y
    c = vcat(center(Z1), center(Z2))
    G = blockdiag(sparse(genmat(Z1)), sparse(genmat(Z2)))
    return Zonotope(c, G)
end
=#

#=
function project!(Z::AbstractZonotope{N}, vars::AbstractVector{Int}) where {N}
    M = projection_matrix(vars, dim(Z), N)
    return linear_map!(M, Z)
end

function project!(X::CartesianProduct{N, <:Interval, <:AbstractZonotope}, vars::AbstractVector{Int}) where {N}
    Z = convert(Zonotope, X)
    return project(Z, vars)
end
=#

#=
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
=#

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

# fallback lazy projection
function _Projection(X::LazySet, vars::NTuple{D, T}) where {D, T<:Integer}
    return LazySets.Projection(X, collect(vars))
end

# ===============================
# Decompositions and partitions
# ===============================

#const Partition{}

# concrete decomposition using a uniform block partition
#using LazySets.Arrays: projection_matrix

#const Partition{N, VT} = AbstractVector{VT} where {VT<:AbstractVector{Int}}

function _decompose(X::LazySet{N},
                    partition::AbstractVector{<:AbstractVector{Int}},
                    set_type::Type{ST}) where {N, ST<:LazySet}
    n = dim(X)
    result = Vector{ST}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        πS = Projection(S, block)
        result[i] = overapproximate(πS, X)
    end
    return CartesianProductArray(result)
end

# =========================
# Overapproximation
# =========================

function _overapproximate(lm::LinearMap{N, <:AbstractZonotope{N}, NM, <:AbstractIntervalMatrix{NM}},
                            ::Type{<:Zonotope}) where {N<:Real, NM}

    Mc, Ms = _split(matrix(lm))
    Z = LazySets.set(lm)
    c = center(Z)
    G = genmat(Z)
    _overapproximate_interval_linear_map(Mc, Ms, c, G)
end

function _overapproximate_interval_linear_map(Mc::AbstractMatrix{N},
                                              Ms::AbstractMatrix{N},
                                              c::AbstractVector{N},
                                              G::AbstractMatrix{N}) where {N}
    n = length(c)
    m = size(G, 2) # number of generators
    c_oa = Mc * c
    Ggens = Mc * G

    Dv = zeros(N, n, n)
    Gt = transpose(G)
    @inbounds for i in 1:n
        Dv[i, i] = abs(c[i])
        for j in 1:m
            Dv[i, i] += abs(Gt[j, i])
        end
    end
    G_oa = hcat(Ggens, Ms * Dv)
    return Zonotope(c_oa, G_oa)
end

function _overapproximate_interval_linear_map(Mc::StaticArray{Tuple{q, n}, T, 2},
                                              Ms::StaticArray{Tuple{q, n}, T, 2},
                                              c::AbstractVector{N},
                                              G::AbstractMatrix{N}) where {N, q, n, T}
    m = size(G, 2) # number of generators
    c_oa = Mc * c
    Ggens = Mc * G

    Dv = zeros(MMatrix{n, n, N})
    c_abs = abs.(c)
    Gt = transpose(G)
    @inbounds for i in 1:n
        Dv[i, i] = c_abs[i]
        for j in 1:m
            Dv[i, i] += abs(Gt[j, i])
        end
    end
    G_oa = hcat(Ggens, Ms * Dv)

    return Zonotope(c_oa, G_oa)
end

function _split_fallback!(A::IntervalMatrix{T}, C, S) where {T}
    m, n = size(A)
    @inbounds for j in 1:n
        for i in 1:m
            itv = A[i, j]
            radius = (sup(itv) - inf(itv)) / T(2)
            C[i, j] = inf(itv) + radius
            S[i, j] = radius
        end
    end
    return C, S
end

function _split(A::IntervalMatrix{T, IT, MT}) where {T, IT, MT<:AbstractMatrix{IT}}
    m, n = size(A)
    C = Matrix{T}(undef, m, n)
    S = Matrix{T}(undef, m, n)
    _split_fallback!(A, C, S)
    return C, S
end

function _split(A::IntervalMatrix{T, IT, MT}) where {T, IT, ST, MT<:StaticArray{ST, IT}}
    m, n = size(A)
    # TODO: use MMatrix and convert to SMatrix afterwards?
    C = Matrix{T}(undef, m, n)
    S = Matrix{T}(undef, m, n)
    _split_fallback!(A, C, S)
    return SMatrix{m, n, T}(C), SMatrix{m, n, T}(S)
end

function _symmetric_interval_hull(x::Interval)
    abs_inf = abs(min(x))
    abs_sup = abs(max(x))
    bound = max(abs_sup, abs_inf)
    return Interval(-bound, bound)
end

# ==================================
# Zonotope order reduction methods
# ==================================

# algorithm selection
#_reduce_order(Z::Zonotope, r, alg::Val{:combastel}) = _reduce_order_COMB03(Z, r)
#_reduce_order(Z::Zonotope, r, alg::Val{:girard}) = _reduce_order_GIR05(Z, r)
_reduce_order(Z, r) = _reduce_order_COMB03(Z, r) # TEMP default > add option to algorithm structs

# Implements zonotope order reduction method from [COMB03]
# We follow the notation from [YS18]
function _reduce_order_COMB03(Z::Zonotope{N, MN}, r::Number) where {N, MN}
    r >= 1 || throw(ArgumentError("the target order should be at least 1, but it is $r"))
    c = Z.center
    G = Z.generators
    n, p = size(G)

    # r is bigger than that order of Z => don't reduce
    (r * n >= p) && return Z

    # compute interval hull of L
    Lred = zeros(N, n, n)
    @inbounds for i in 1:n
        for j in 1:p
            Lred[i, i] += abs(G[i, j])
        end
    end

    if isone(r)
        return Zonotope(c, Lred)
    end

    # split Z into K and L
    # sort generators by decreasing 2-norm
    weights = zeros(N, p)
    @inbounds for i in 1:p
        weights[i] = norm(view(G, :, i), 2)
    end
    indices = Vector{Int}(undef, p)
    sortperm!(indices, weights, rev=true, initialized=false)

    m = floor(Int, n * (r - 1))
    Gred = hcat(G[:, indices[1:m]], Lred)
    return Zonotope(c, Gred)
end

# Implements zonotope order reduction method from [GIR05]
# We follow the notation from [YS18]
function _reduce_order_GIR05(Z::Zonotope{N, MN}, r::Number) where {N, MN}
    r >= 1 || throw(ArgumentError("the target order should be at least 1, but it is $r"))
    c = Z.center
    G = Z.generators
    n, p = size(G)

    # r is bigger than that order of Z => don't reduce
    (r * n >= p) && return Z

    # compute interval hull of L
    Lred = zeros(N, n, n)
    @inbounds for i in 1:n
        for j in 1:p
            Lred[i, i] += abs(G[i, j])
        end
    end

    if isone(r)
        return Zonotope(c, Lred)
    end

    # split Z into K and L
    # sort generators by decreasing 2-norm
    weights = zeros(N, p)
    @inbounds for i in 1:p
        v = view(G, :, i)
        weights[i] = norm(v, 1) - norm(v, Inf)
    end
    indices = Vector{Int}(undef, p)
    sortperm!(indices, weights, rev=true, initialized=false)

    m = floor(Int, n * (r - 1))
    Gred = hcat(G[:, indices[1:m]], Lred)
    return Zonotope(c, Gred)
end

# ====================================
# Intersection
# ====================================

using LazySets: _geq, _leq

# fallback
function _is_intersection_empty(X, Y)
    return LazySets.is_intersection_empty(X, Y)
end

function _is_intersection_empty(I1::Interval{N}, I2::Interval{N}) where {N<:Real}
    return !_leq(min(I2), max(I1))  || !_leq(min(I1), max(I2))
end

# H : {x : ax <= b}, one-dimensional with a != 0
function _is_intersection_empty(X::Interval{N}, H::HalfSpace{N}) where {N<:Real}
    a = H.a[1]
    b = H.b
    if a > zero(N)
        return !_leq(min(X), b/a)
    else
        return !_geq(max(X), b/a)
    end
end
# symmetric case
_is_intersection_empty(H::HalfSpace, X::Interval) = _is_intersection_empty(X, H)

# H : {x : ax = b}, one-dimensional with a != 0
function _is_intersection_empty(X::Interval{N}, H::Hyperplane{N}) where {N<:Real}
    q = H.b / H.a[1]
    return !_geq(q, min(X)) || !_leq(q, max(X))
end
# symmetric case
_is_intersection_empty(H::Hyperplane, X::Interval) = _is_intersection_empty(X, H)
