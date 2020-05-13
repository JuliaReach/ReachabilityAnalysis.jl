# =========================
# Conversion
# =========================

@inline nodim_msg() = throw(ArgumentError("to use the `static` option you should pass " *
                                          "the system's dimension argument `dim=...`"))

# no-op
_reconvert(Ω0::Zonotope{N, Vector{N}, Matrix{N}}, static::Val{false}, dim, ngens) where {N} = Ω0
_reconvert(Ω0::Zonotope{N, <:SVector, <:SMatrix}, static::Val{true}, dim) where {N} = Ω0

# convert any zonotope to be represented wih regular arrays
_reconvert(Ω0::Zonotope, static::Val{false}, dim, ngens) = Zonotope(Vector(Ω0.center), Matrix(Ω0.generators))

# convert any zonotope to be represented with static arrays
function _reconvert(Ω0::Zonotope{N, VN, MN}, static::Val{true}, dim::Val{n}, ngens::Val{p}) where {N, VN, MN, n, p}
    G = Ω0.generators
    m = size(G, 2)
    c = SVector{n, N}(Ω0.center)

    if m == p
        return Zonotope(c, SMatrix{n, p}(G))

    elseif m < p
        # extend with zeros
        Gext = hcat(SMatrix{n, m}(G), zeros(MMatrix{n, p-m, N}))
        return Zonotope(c, Gext)

    else
        throw(ArgumentError("can't reconvert a zontope with $m generators to a " *
                            "zonotope with $p generators; you should reduce it first"))
    end
end

# no-op
_reconvert(Ω0::Hyperrectangle{N, Vector{N}, Vector{N}}, static::Val{false}, dim::Missing) where {N} = Ω0
_reconvert(Ω0::Hyperrectangle{N, Vector{N}, Vector{N}}, static::Val{true}, dim::Missing) where {N} = nodim_msg()
_reconvert(Ω0::Hyperrectangle{N, <:SVector, <:SVector}, static::Val{true}, dim) where {N} = Ω0

# convert any Hyperrectangle to be represented wih regular arrays
function _reconvert(Ω0::Hyperrectangle, static::Val{false}, dim::Missing)
    Ω0 = Hyperrectangle(Vector(Ω0.center), Matrix(Ω0.radius), check_bounds=false)
end

# convert any Hyperrectangle to be represented with static arrays
function _reconvert(Ω0::Hyperrectangle{N, VNC, VNR}, static::Val{true}, dim::Val{n}) where {N, VNC, VNR, n}
    #n = length(Ω0.center) # dimension
    Ω0 = Hyperrectangle(SVector{n, N}(Ω0.center), SVector{n, N}(Ω0.radius), check_bounds=false)
end

# no-op
_reconvert(Φ::Matrix{N}, static::Val{false}, dim) where {N} = Φ
_reconvert(Φ::IntervalMatrix{N}, static::Val{false}, dim) where {N} = Φ
_reconvert(Φ::AbstractMatrix, static::Val{false}, dim) = Matrix(Φ)
_reconvert(Φ::SMatrix, static::Val{true}, dim) = Φ
_reconvert(Φ::AbstractMatrix, static::Val{true}, dim::Missing) = nodim_msg()

function _reconvert(Φ::AbstractMatrix{N}, static::Val{true}, dim::Val{n}) where {N, n}
    #n = size(Φ, 1)
    Φ = SMatrix{n, n, N, n*n}(Φ)
end

function _reconvert(Φ::IntervalMatrix{N, IN, Matrix{IN}}, static::Val{true}, dim::Val{n}) where {N, IN, n}
    #n = size(Φ, 1)
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

function _convert_or_overapproximate(X::LazySet, T::Type{<:AbstractPolytope})
    return _convert_or_overapproximate(T, X)
end

# =========================
# In-place ops
# =========================

# Extension of some common LazySets operations, some of them in-place

# TODO: remove
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

# TEMP
@inline _linear_map(M, X) = LazySets.linear_map(M, X)

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
                    blocks, # ::AbstractVector{<:AbstractVector{Int}}
                    set_type::Type{ST}) where {N, ST<:LazySet}
    n = dim(X)
    result = Vector{ST}(undef, length(blocks))

    @inbounds for (i, bi) in enumerate(blocks)
        πX = _Projection(X, bi)
        result[i] = overapproximate(πX, ST)
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

    dvec = zeros(N, n)
    @inbounds for i in 1:n
        dvec[i] = abs(c[i])
        for j in 1:m
            dvec[i] += abs(G[i, j])
        end
    end
    DV = zeros(N, n, n)
    α = Ms * dvec
    @inbounds for i in 1:n
        DV[i, i] = α[i]
    end
    G_oa = hcat(Ggens, DV)
    return Zonotope(c_oa, G_oa)
end

function _overapproximate_interval_linear_map(Mc::SMatrix{n, n, N, LM},
                                              Ms::SMatrix{n, n, N, LM},
                                              c::SVector{n, N},
                                              G::SMatrix{n, m, N, LG}) where {n, N, LM, m, LG}
    c_oa = Mc * c
    Ggens = Mc * G

    dvec = zeros(N, n)
    @inbounds for i in 1:n
        dvec[i] = abs(c[i])
        for j in 1:m
            dvec[i] += abs(G[i, j])
        end
    end
    DV = zeros(MMatrix{n, n, N}) # NOTE: sole difference with regular arrays, may refactor
    α = Ms * dvec
    @inbounds for i in 1:n
        DV[i, i] = α[i]
    end
    G_oa = hcat(Ggens, DV)
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

# attempts type-stability
function _symmetric_interval_hull(x::Interval)
    abs_inf = abs(min(x))
    abs_sup = abs(max(x))
    bound = max(abs_sup, abs_inf)
    return Interval(-bound, bound)
end

# attempts type-stability
function _symmetric_interval_hull(S::LazySet{N}) where {N<:Real}
    # fallback returns a hyperrectangular set
    (c, r) = LazySets.Approximations.box_approximation_helper(S)
    #if r[1] < 0
    #    return EmptySet{N}(dim(S))
    #end
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

# attemp type stability
function _overapproximate(S::LazySet{N}, ::Type{<:Hyperrectangle}) where {N<:Real}
    c, r = LazySets.Approximations.box_approximation_helper(S)
    #if r[1] < 0
    #    return EmptySet{N}(dim(S))
    #end
    return Hyperrectangle(c, r)
end

# TEMP
function LazySets.Approximations.box_approximation(x::IntervalArithmetic.Interval)
    return convert(Hyperrectangle, Interval(x))
end

# TEMP
function LazySets.Approximations.box_approximation(x::IntervalArithmetic.IntervalBox)
    return convert(Hyperrectangle, x)
end

# concrete set complement for polyhedral sets
complement(X::LazySet) = UnionSetArray(constraints_list(Complement(X)))

# list of constraints of the set complement of a polyhedral set
function LazySets.constraints_list(CX::Complement{N, ST}) where {N, ST<:AbstractPolyhedron{N}}
    clist = constraints_list(CX.X)
    out = similar(clist)
    for (i, ci) in enumerate(clist)
        out[i] = LinearConstraint(-ci.a, -ci.b)
    end
    return out
end

function _overapproximate(S::UnionSetArray{N}, ::Type{<:Hyperrectangle}) where {N}
    c, r = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}(dim(S))
    end
    return Hyperrectangle(c, r)
end

@inline function box_approximation_helper(S::UnionSetArray{N}) where {N}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
        if r[i] < 0
            # contradicting bounds => set is empty
            # terminate with first radius entry being negative
            r[1] = r[i]
            break
        end
    end
    return c, r
end

# ==================================
# Zonotope splitting methods
# ==================================

function _split(Z::Zonotope{N, Vector{N}, Matrix{N}}, j::Int) where {N}
    c, G = Z.center, Z.generators
    n, p = size(G)
    @assert 1 <= j <= p "cannot split a zonotope with $p generators along index $j"

    c₁ = Vector{N}(undef, n)
    c₂ = Vector{N}(undef, n)

    G₁ = copy(G)
    G₂ = copy(G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
        G₁[i, j] = α
    end

    Z₁ = Zonotope(c₁, G₁)
    Z₂ = Zonotope(c₂, G₂)
    return Z₁, Z₂
end

function _split(Z::Zonotope{N, SVector{n, N}, <:SMatrix{n, p, N}}, j::Int) where {N, n, p}
    @assert 1 <= j <= p "cannot split a zonotope with $p generators along index $j"
    c, G = Z.center, Z.generators

    c₁ = MVector{n, N}(undef)
    c₂ = MVector{n, N}(undef)

    G₁ = MMatrix{n, p}(G)
    G₂ = MMatrix{n, p}(G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
        G₁[i, j] = α
    end

    Z₁ = Zonotope(SVector{n}(c₁), SMatrix{n, p}(G₁))
    Z₂ = Zonotope(SVector{n}(c₂), SMatrix{n, p}(G₂))
    return Z₁, Z₂
end

# ==================================
# Zonotope order reduction methods
# ==================================

abstract type AbstractReductionMethod end

# These methods split the given zonotope Z into two zonotopes, K and L, where
# K contains the most "representative" generators and L contains the generators
# that are reduced, Lred using a box overapproximation
struct GIR05 <: AbstractReductionMethod end
struct COMB03 <: AbstractReductionMethod end

const _COMB03 = COMB03()
const _GIR05 = GIR05()

# algorithm selection
_reduce_order(Z::Zonotope, r::Number) = _reduce_order_GIR05(Z, r) # default
_reduce_order(Z::Zonotope, r::Number, ::GIR05) = _reduce_order_GIR05(Z, r)
_reduce_order(Z::Zonotope, r::Number, ::COMB03) = _reduce_order_COMB03(Z, r)

# return the indices of the generators in G (= columns) sorted according to the COMB03 method
# the generator index with highest score goes first
function _weighted_gens!(indices, G::AbstractMatrix{N}, ::COMB03) where {N}
    p = size(G, 2)
    weights = Vector{N}(undef, p)
    @inbounds for j in 1:p
        v = view(G, :, j)
        weights[j] = norm(v, 2)
    end
    sortperm!(indices, weights, rev=true, initialized=false)
    return indices
end

# return the indices of the generators in G (= columns) sorted according to the GIR05 method
# the generator index with highest score goes first
function _weighted_gens!(indices, G::AbstractMatrix{N}, ::GIR05) where {N}
    n, p = size(G)
    weights = Vector{N}(undef, p)
    @inbounds for j in 1:p
        aux_norm_1 = zero(N)
        aux_norm_inf = zero(N)
        for i in 1:n
            abs_Gij = abs(G[i, j])
            aux_norm_1 += abs_Gij
            if aux_norm_inf < abs_Gij
                aux_norm_inf = abs_Gij
            end
        end
        weights[j] = aux_norm_1 - aux_norm_inf
    end
    sortperm!(indices, weights, rev=true, initialized=false)
    return indices
end

# compute interval hull of the generators of G (= columns) corresponding to `indices`
function _interval_hull(G::AbstractMatrix{N}, indices) where {N}
    n, p = size(G)
    Lred = zeros(N, n, n)
    @inbounds for i in 1:n
        for j in indices
            Lred[i, i] += abs(G[i, j])
        end
    end
    return Lred
end

# implementation for static arrays
function _interval_hull(G::SMatrix{n, p, N, L}, indices) where {n, p, N, L}
    Lred = zeros(MMatrix{n, n, N})
    @inbounds for i in 1:n
        for j in indices
            Lred[i, i] += abs(G[i, j])
        end
    end
    return SMatrix{n, n}(Lred)
end

# given an n x p matrix G and a vector of m integer indices with m <= p,
# concatenate the columns of G given by `indices` with the matrix Lred
function _hcat_KLred(G::AbstractMatrix, indices, Lred::AbstractMatrix)
    K = view(G, :, indices)
    return hcat(K, Lred)
end

# implementation for static arrays
function _hcat_KLred(G::SMatrix{n, p, N, L1}, indices, Lred::SMatrix{n, n, N, L2}) where {n, p, N, L1, L2}
    m = length(indices)
    K = SMatrix{n, m}(view(G, :, indices))
    return hcat(K, Lred)
end

# Implements zonotope order reduction method from [COMB03]
# We follow the notation from [YS18]
function _reduce_order_COMB03(Z::Zonotope{N}, r::Number) where {N}
    r >= 1 || throw(ArgumentError("the target order should be at least 1, but it is $r"))
    c = Z.center
    G = Z.generators
    n, p = size(G)

    # r is bigger than the order of Z => don't reduce
    (r * n >= p) && return Z

    # this algorithm sort generators by decreasing 2-norm
    indices = Vector{Int}(undef, p)
    _weighted_gens!(indices, G, _COMB03)

    # the first m generators have greatest weight
    m = floor(Int, n * (r - 1))

    # compute interval hull of L
    Lred = _interval_hull(G, view(indices, (m+1):p))

    isone(r) && return Zonotope(c, Lred)

    Gred = _hcat_KLred(G, view(indices, 1:m), Lred)
    return Zonotope(c, Gred)
end

# Implements zonotope order reduction method from [GIR05]
# We follow the notation from [YS18]
function _reduce_order_GIR05(Z::Zonotope{N}, r::Number) where {N}
    r >= 1 || throw(ArgumentError("the target order should be at least 1, but it is $r"))
    c = Z.center
    G = Z.generators
    n, p = size(G)

    # r is bigger than the order of Z => don't reduce
    (r * n >= p) && return Z

    # this algorithm sorts generators by ||⋅||₁ - ||⋅||∞ difference
    indices = Vector{Int}(undef, p)
    _weighted_gens!(indices, G, _GIR05)

    # the first m generators have greatest weight
    m = floor(Int, n * (r - 1))

    # compute interval hull of L
    Lred = _interval_hull(G, view(indices, (m+1):p))

    isone(r) && return Zonotope(c, Lred)

    Gred = _hcat_KLred(G, view(indices, 1:m), Lred)
    return Zonotope(c, Gred)
end

# ====================================
# Disjointness checks
# ====================================

abstract type AbstractDisjointnessMethod end

struct ZonotopeEnclosure <: AbstractDisjointnessMethod
# we overapproximate the reach-set with a zonotope, then make the disjointness check
end

struct BoxEnclosure <: AbstractDisjointnessMethod
# we overapproximate the reach-set with a hyperrectangle, then make the disjointness check
end

struct Exact <: AbstractDisjointnessMethod
# we pass the sets to is_intersection_empty without pre-processing
end

# fallback
_is_intersection_empty(X::LazySet, Y::LazySet) = LazySets.is_intersection_empty(X, Y)

# symmetric case
_is_intersection_empty(X::LazySet, R::AbstractReachSet, method=Exact()) = _is_intersection_empty(R, X, method)

# --------------------------------------------------------------------
# Methods to evaluate disjointness between a reach-set and a lazy set
# --------------------------------------------------------------------

function _is_intersection_empty(X::AbstractReachSet, Y::LazySet, ::ZonotopeEnclosure)
    Z = overapproximate(X, Zonotope)
    return _is_intersection_empty(set(Z), Y)
end

function _is_intersection_empty(X::AbstractReachSet, Y::LazySet, ::BoxEnclosure)
    H = overapproximate(X, Hyperrectangle)
    return _is_intersection_empty(set(H), Y)
end

# -----------------------------------------------
# Disjointness checks between specific set types
# -----------------------------------------------

using LazySets: _geq, _leq

function _is_intersection_empty(I1::Interval{N}, I2::Interval{N}) where {N<:Real}
    return !_leq(min(I2), max(I1)) || !_leq(min(I1), max(I2))
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

# ==================================
# Concrete intersection
# ==================================

# ---------------------------------------------------------
# Methods to evaluate intersection between a pair of sets
# ---------------------------------------------------------

abstract type AbstractIntersectionMethod end

struct BoxEnclosure <: AbstractIntersectionMethod
#
end

struct Exact <: AbstractIntersectionMethod
#
end
