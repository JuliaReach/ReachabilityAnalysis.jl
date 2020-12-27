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
function _reconvert(Ω0::Zonotope{N, VN, MN}, static::Val{true}, dim::Missing, ngens::Missing) where {N, VN, MN}
    n, m = size(Ω0.generators)
    _reconvert(Ω0, static, Val(n), Val(m))
end

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
_reconvert(Φ::AbstractMatrix, static::Val{true}, dim::Missing) = _reconvert(Φ, static, Val(size(Φ, 1)))

function _reconvert(Φ::AbstractMatrix{N}, static::Val{true}, dim::Val{n}) where {N, n}
    @assert n == size(Φ, 1)
    Φ = SMatrix{n, n, N, n*n}(Φ)
end

function _reconvert(Φ::IntervalMatrix{N, IN, Matrix{IN}}, static::Val{true}, dim::Val{n}) where {N, IN, n}
    @assert n == size(Φ, 1)
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

function _overapproximate(X::Hyperrectangle, T::Type{HPolytope{N, VT}}) where {N, VT}
    # TODO create overapproximation using VT directly
    Y = overapproximate(X, BoxDirections(dim(X)))
    return convert(T, Y)
end

Base.convert(::Type{HPolytope{N, VT}}, P::HPolytope{N, VT}) where {N, VT} = P
Base.convert(::Type{HPolytope{N, VT}}, P) where {N, VT} = HPolytope([HalfSpace(VT(c.a), c.b) for c in constraints_list(P)])

function Base.convert(::Type{Hyperrectangle{N, Vector{N}, Vector{N}}},
                      H::Hyperrectangle{N, SVector{L, N}, SVector{L, N}}) where {N, L}
    return Hyperrectangle(Vector(H.center), Vector(H.radius))
end

function Base.convert(::Type{Singleton},
                      cp::CartesianProduct{N, S1, S2}) where{N, S1<:Singleton{N}, S2<:Singleton{N}}
    x = element(cp.X)
    y = element(cp.Y)
    return Singleton(vcat(x, y))
end

# =========================
# Projection
# =========================

# extend LazySets concrete projection for other arg fomats
LazySets.project(X::LazySet, vars::NTuple{D, <:Integer}) where {D} = project(X, collect(vars))
LazySets.project(X::LazySet; vars) where {D} = project(X, collect(vars))

# extend LazySets lazy projection for other arg fomats
LazySets.Projection(X::LazySet, vars::NTuple{D, <:Integer}) where {D} = Projection(X, collect(vars))
LazySets.Projection(X::LazySet; vars) where {D} = Projection(X, collect(vars))

# fallback concrete projection
function _project(X::LazySet, vars::NTuple{D, <:Integer}) where {D}
    return project(X, collect(vars))
end

_project(X::LazySet, vars::AbstractVector) = _project(X, Tuple(vars))

# more efficient concrete projection for zonotopic sets (see LazySets#2013)
function _project(Z::AbstractZonotope{N}, vars::NTuple{D, <:Integer}) where {N, D}
    M = projection_matrix(collect(vars), dim(Z), N)
    return linear_map(M, Z)
end

# cartesian product of an interval and a hyperrectangular set
# TODO refactor to LazySets
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
# TODO refactor after LazySets#2190
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

# fallback lazy projection
function _Projection(X::LazySet, vars::NTuple{D, <:Integer}) where {D}
    return LazySets.Projection(X, collect(vars))
end

_Projection(X::LazySet, vars::AbstractVector) =  _Projection(X, Tuple(vars))

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
    q = size(Mc, 1)
    α = Ms * dvec # vector of length q
    αnz = findall(!iszero, α)
    DV = zeros(N, q, length(αnz))
    @inbounds for (j, idx) in enumerate(αnz)
        DV[j, idx] = α[idx]
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
    q = size(Mc, 1)
    α = Ms * dvec # vector of length q
    αnz = findall(!iszero, α)
    DV = zeros(MMatrix{q, q, N}) # NOTE: sole difference with regular arrays, may refactor
    @inbounds for (j, idx) in enumerate(αnz)
        DV[j, idx] = α[idx]
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

# fix return type to hyperrectangle for type stability
function _overapproximate(S::LazySet{N}, ::Type{<:Hyperrectangle}) where {N<:Real}
    c, r = LazySets.Approximations.box_approximation_helper(S)
    #if r[1] < 0
    #    return EmptySet{N}(dim(S))
    #end
    return Hyperrectangle(c, r)
end

function LazySets.Approximations.box_approximation(x::IntervalArithmetic.Interval)
    return convert(Hyperrectangle, Interval(x))
end

function LazySets.Approximations.box_approximation(x::IntervalArithmetic.IntervalBox)
    return convert(Hyperrectangle, x)
end

# concrete set complement for polyhedral sets
# See LazySets#2381
complement(H::HalfSpace) = HalfSpace(-H.a, -H.b)
complement(X::LazySet) = UnionSetArray(constraints_list(Complement(X)))

# list of constraints of the set complement of a polyhedral set
# TODO refactor to LazySets
function LazySets.constraints_list(CX::Complement{N, ST}) where {N, ST<:AbstractPolyhedron{N}}
    clist = constraints_list(CX.X)
    out = similar(clist)
    for (i, ci) in enumerate(clist)
        out[i] = LinearConstraint(-ci.a, -ci.b)
    end
    return out
end

LazySets.box_approximation(S::UnionSetArray) = overapproximate(S, Hyperrectangle)

function LazySets.overapproximate(S::UnionSetArray{N}, ::Type{<:Hyperrectangle}) where {N}
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

# overapproximate a hyperrectangular set with a polytope
_overapproximate(H::AbstractHyperrectangle, T::Type{<:HPolytope}) = _overapproximate_hyperrectangle(H, T)
_overapproximate(H::Hyperrectangle, T::Type{<:HPolytope}) = _overapproximate_hyperrectangle(H, T)
function _overapproximate_hyperrectangle(H, ::Type{<:HPolytope})
    P = overapproximate(H, Hyperrectangle)
    HPolytope([HalfSpace(Vector(c.a), c.b) for c in constraints_list(H)])
end

"""
    relative_error(x, x_ref)

Compute the relative error between interval `x` and a reference interval `xref`.

### Input

- `x`    -- interval
- `xref` -- reference interval

### Output

An interval representing the relative error (in percentage) of `x` with respect to
the reference interval `xref`.

### Algorithm

If ``x = [x_L, x_H]``` and ``xref = [xref_L, xref_H]``, the output is the interval
``y = 100 * [y_L, y_H]`` computed as ``y_L = -(x_L - xref_L) / den`` and
``y_H = (x_H - xref_H) / den``, where ``den = xref_H - xref_L``.

This function measures the relative error between an interval `x` and a reference
interval `x_ref` accounting for it the lower and the upper range bounds separately
(see  Eq. (20) in [1]).

### References

- [1] Althoff, Matthias, Dmitry Grebenyuk, and Niklas Kochdumper.
      "Implementation of Taylor models in CORA 2018."
      Proc. of the 5th International Workshop on Applied Verification for Continuous
      and Hybrid Systems. 2018. [pdf](https://easychair.org/publications/paper/9Tz3)
"""
function relative_error(x, x_ref)
    x_low, x_high = inf(x), sup(x)
    x_ref_low, x_ref_high = inf(x_ref), sup(x_ref)
    denom = x_ref_high - x_ref_low
    rel_low = -(x_low - x_ref_low) / denom
    rel_high = (x_high - x_ref_high) / denom
    return 100 * IntervalArithmetic.Interval(rel_low, rel_high)
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

# we overapproximate the reach-set with a zonotope, then make the disjointness check
struct ZonotopeEnclosure <: AbstractDisjointnessMethod end

# we overapproximate the reach-set with a hyperrectangle, then make the disjointness check
struct BoxEnclosure <: AbstractDisjointnessMethod end

# we pass the sets to the disjointness check algorithm without pre-processing
struct NoEnclosure <: AbstractDisjointnessMethod end

# this is a dummy disjointness check which returns "false" irrespective of the value of its arguments
struct Dummy <: AbstractDisjointnessMethod end

# LazySets fallback
_is_intersection_empty(X::LazySet, Y::LazySet, ::NoEnclosure) = is_intersection_empty(X, Y)

# enclose the first argument with a zonotope
function _is_intersection_empty(X::LazySet, Y::LazySet, ::ZonotopeEnclosure)
    Z = overapproximate(X, Zonotope)
    return is_intersection_empty(Z, Y)
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, method::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return _is_intersection_empty(set(Z), Y)
end

@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, method::BoxEnclosure)
    H = overapproximate(R, Hyperrectangle)
    return _is_intersection_empty(set(H), Y)
end

# in this method we assume that the intersection is non-empty
@commutative function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, method::Dummy)
    return false
end

# TODO refactor?
function overapproximate(R::TaylorModelReachSet{N}, ::BoxEnclosure) where {N}
    return overapproximate(R, Hyperrectangle)
end
# TODO refactor?
function overapproximate(R::TaylorModelReachSet{N}, ::ZonotopeEnclosure) where {N}
    return overapproximate(R, Zonotope)
end

# -----------------------------------------------
# Disjointness checks between specific set types
# -----------------------------------------------

using LazySets: _geq, _leq

# H : {x : ax <= b}, one-dimensional with a != 0
@commutative function _is_intersection_empty(X::Interval, H::HalfSpace, method::NoEnclosure)
    a = H.a[1]
    b = H.b
    N = promote_type(eltype(X), eltype(H))
    if a > zero(N)
        return !_leq(min(X), b/a)
    else
        return !_geq(max(X), b/a)
    end
end

# H : {x : ax = b}, one-dimensional with a != 0
@commutative function _is_intersection_empty(X::Interval, H::Hyperplane, method::NoEnclosure)
    q = H.b / H.a[1]
    return !_geq(q, min(X)) || !_leq(q, max(X))
end

# if X is polyhedral and Y is the set union of half-spaces, X ∩ Y is empty iff
# X ∩ Hi is empty for each half-space Hi in Y
# NOTE the algorithm below solves an LP for each X ∩ Hi; however, we can proceed
# more efficintly using support functions
# see LazySets.is_intersection_empty_helper_halfspace
function LazySets.is_intersection_empty(X::AbstractPolytope{N},
                                        Y::UnionSetArray{N, HT}
                                        ) where {N<:Real, VN<:AbstractVector{N}, HT<:HalfSpace{N, VN}}
    clist_X = constraints_list(X)
    for ci in Y.array
        # using support functions
        #!(-ρ(-hs.a, X) > hs.b) && return false # TODO use robust check

        # using LP
        Y_ci = vcat(clist_X, ci)
        remove_redundant_constraints!(Y_ci) && return false
    end
    return true
end

# ==================================
# Concrete intersection
# ==================================

# -------------------------------------------------------------
# Methods to compute the intersection between two or more sets
# ------------------------------------------------------------

abstract type AbstractIntersectionMethod end

#= TODO annotate normal vector types?
struct HRepIntersection{SX, SY} <: AbstractIntersectionMethod end
setrep(int_method::HRepIntersection{SX<:AbstractPolytope}, SY<:AbstractPolyhedron}) = SX
=#

# evaluate X ∩ Y exactly using the constraint representation of X and Y
# evaluate X₁ ∩ ⋯ ∩ Xₖ using the constraint representation of each Xᵢ
struct HRepIntersection <: AbstractIntersectionMethod
#
end

setrep(::HRepIntersection) = HPolytope{Float64, Vector{Float64}}

# "fallback" implementation that uses LazySets intersection(X, Y)
struct FallbackIntersection <: AbstractIntersectionMethod
#
end

struct BoxIntersection <: AbstractIntersectionMethod
#
end

setrep(::BoxIntersection) = Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}

# evaluate X ∩ Y approximately using support function evaluations, with the
# property that the support function of X ∩ Y along direction d is not greater
# ρ(d, X ∩ Y) <= min(ρ(d, X), ρ(d, Y))
# by the same token, compute X₁ ∩ ⋯ ∩ Xₖ approximately using the same property
# if the template is provided, we have TN<:AbstractDirections{N, VN}
# otherwise, the constraints of X and Y are used and TN is Missing
struct TemplateHullIntersection{N, VN, TN} <: AbstractIntersectionMethod
    dirs::TN
end

# constructor with template directions provided
function TemplateHullIntersection(dirs::TN) where {N, VN, TN<:AbstractDirections{N, VN}}
    TemplateHullIntersection{N, VN, TN}(dirs)
end

# constructor without template directions => directions are missing until evaluated
function TemplateHullIntersection{N, VN}() where {N, VN<:AbstractVector{N}}
    TemplateHullIntersection{N, VN, Missing}(missing)
end
TemplateHullIntersection() = TemplateHullIntersection{Float64, Vector{Float64}}()

setrep(::TemplateHullIntersection{N, VN}) where {N, VN} = HPolytope{N, VN}
setrep(::TemplateHullIntersection{N, SEV}) where {N, SEV<:SingleEntryVector{N}} = Union{HPolytope{N, SEV}, HPolytope{N, Vector{N}}}
setrep(::TemplateHullIntersection{N, SP}) where {N, SP<:SparseVector{N}} = Union{HPolytope{N, SP}, HPolytope{N, Vector{N}}}

_intersection(X::LazySet, Y::LazySet, ::FallbackIntersection) = intersection(X, Y)
_intersection(X::LazySet, Y::LazySet) = _intersection(X, Y, FallbackIntersection())

# converts the normal vector of a list of half-spaces to be a Vector
const VECH{N, VT} = Vector{HalfSpace{N, VT}}
_to_vec(c::HalfSpace{N, Vector{N}}) where {N} = c
_to_vec(c::HalfSpace{N, VT}) where {N, VT<:AbstractVector{N}} = HalfSpace(Vector(c.a), c.b)
_to_vec(x::VECH{N, Vector{N}}) where {N} = x
_to_vec(x::VECH{N, VT}) where {N, VT<:AbstractVector{N}} = broadcast(_to_vec, x)

# concatenates lists of half-spaces such that the normal vectors of the final list
# are all Vector
_vcat(args::VECH...) = vcat(map(_to_vec, args)...)

# TODO use LazySets intersection
function _intersection(X::AbstractPolyhedron, Y::AbstractPolyhedron, ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    out = _vcat(clist_X, clist_Y)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

function _intersection(X::AbstractPolyhedron, Y::AbstractPolyhedron, Z::AbstractPolyhedron, ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    clist_Z = constraints_list(Z)
    out = _vcat(clist_X, clist_Y, clist_Z)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

function _intersection(X::LazySet, Y::AbstractPolyhedron, Z::AbstractPolyhedron, W::AbstractPolyhedron, ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    clist_Z = constraints_list(Z)
    clist_W = constraints_list(W)
    out = _vcat(clist_X, clist_Y, clist_Z, clist_W)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

# if the template directions is missing => use the constraints of both X and Y
# TODO remove redundant constraints?
function _intersection(X::LazySet, Y::LazySet, method::TemplateHullIntersection{N, VN, Missing}) where {N, VN}
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)

    out = Vector{HalfSpace{N, VN}}()
    sizehint!(out, length(clist_X) + length(clist_Y))

    @inbounds for (i, c) in enumerate(clist_X)
        d = convert(VN, c.a) # normal vector
        b = min(c.b, ρ(d, Y)) # we know that ρ(d, X) = b
        push!(out, HalfSpace(d, b))
    end

    @inbounds for (i, c) in enumerate(clist_Y)
        d = convert(VN, c.a) # normal vector
        b = min(ρ(d, X), c.b) # we know that ρ(d, Y) = b
        push!(out, HalfSpace(d, b))
    end

    return HPolytope(out)
end

# use ρ(d, X∩Y) ≤ min(ρ(d, X), ρ(d, Y)) for each d in the template
function _intersection(X::LazySet, Y::LazySet, method::TemplateHullIntersection{N, VN, TN}) where {N, VN, TN<:AbstractDirections{N, VN}}
    dirs = method.dirs
    out = Vector{HalfSpace{N, VN}}(undef, length(dirs))
    @inbounds for (i, d) in enumerate(dirs)
        d = convert(VN, d)
        b = min(ρ(d, X), ρ(d, Y))
        out[i] = HalfSpace(d, b)
    end
    return HPolytope(out)
end

# =====================================
# Methods for checking inclusion
# =====================================

# given two reach-sets A, B check whether f(A) ⊆ g(B) where
# A ⊆ f(A) and B ⊆ g(B)
abstract type AbstractInclusionMethod end

# no pre-processing of the sets
struct FallbackInclusion <: AbstractInclusionMethod end
