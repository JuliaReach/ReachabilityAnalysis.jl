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
_Projection(X::LazySet, vars::AbstractVector) = LazySets.Projection(X, vars)

# using a vars tuple
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

# TODO pass backend
function LazySets.vertices_list(X::ConvexHullArray{N, PT}) where {N, PT<:AbstractPolytope{N}}
    vlist = [vertices_list(Xi) for Xi in array(X)]
    return convex_hull(vcat(vlist...))
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

# we overapproximate the reach-set with a zonotope, then make the disjointness check
struct ZonotopeEnclosure <: AbstractDisjointnessMethod end

# we overapproximate the reach-set with a hyperrectangle, then make the disjointness check
struct BoxEnclosure <: AbstractDisjointnessMethod end

# we pass the sets to the disjointness check algorithm without pre-processing
struct NoEnclosure <: AbstractDisjointnessMethod end

# fallbacks
_is_intersection_empty(X::LazySet, Y::LazySet; kwargs...) = LazySets.is_intersection_empty(X, Y; kwargs...)
_is_intersection_empty(R::AbstractReachSet, Y::LazySet) = _is_intersection_empty(R, Y, NoEnclosure())

# symmetric case
_is_intersection_empty(X::LazySet, R::AbstractReachSet, method=NoEnclosure()) = _is_intersection_empty(R, X, method)

# this is a dummy disjointness check which returns "false" irrespective of the value of its arguments
struct Dummy <: AbstractDisjointnessMethod end

# --------------------------------------------------------------------
# Methods to evaluate disjointness between a reach-set and a lazy set
# --------------------------------------------------------------------

function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, ::NoEnclosure)
    return _is_intersection_empty(set(R), Y)
end

function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, ::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return _is_intersection_empty(set(Z), Y)
end

function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, ::BoxEnclosure)
    H = overapproximate(R, Hyperrectangle)
    return _is_intersection_empty(set(H), Y)
end

# in this method we assume that the intersection is non-empty
function _is_intersection_empty(R::AbstractReachSet, Y::LazySet, ::Dummy)
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

# propagate methods from reach-set to sets
# TODO always return ReachSets; extend to AbstractLazyReachSet; intersect time spans
_intersection(R::AbstractReachSet, X::LazySet, method::AbstractIntersectionMethod) = _intersection(set(R), X, method)
_intersection(X::LazySet, R::AbstractReachSet, method::AbstractIntersectionMethod) = _intersection(X, set(R), method)

# TODO intersect time spans?
#_intersection(R::AbstractReachSet, S::AbstractReachSet, method::AbstractIntersectionMethod) = _intersection(set(R), set(S), method)

_intersection(X::LazySet, Y::LazySet, ::FallbackIntersection) = intersection(X, Y)
_intersection(X::LazySet, Y::LazySet) = _intersection(X, Y, FallbackIntersection())

# for TMs, we overapproximate with a zonotope
function _intersection(R::TaylorModelReachSet, Y::LazySet)
    _intersection(R, Y, FallbackIntersection())
end

function _intersection(R::TaylorModelReachSet, Y::LazySet, ::FallbackIntersection)
    X = set(overapproximate(R, Zonotope))
    intersection(X, Y)
end

# TODO overapprox Y with a box if it's bounded?
function _intersection(R::TaylorModelReachSet, Y::LazySet, ::BoxIntersection)
    X = set(overapproximate(R, Hyperrectangle))
    out = intersection(X, Y)
    return isempty(out) ? EmptySet(dim(out)) : overapproximate(out, Hyperrectangle)
end

# TODO refactor? See LazySets#2158
function _intersection(R::TaylorModelReachSet, Y::UnionSetArray{N, HT}, ::BoxIntersection) where {N, HT<:HalfSpace{N}}
    X = set(overapproximate(R, Hyperrectangle))

    # find first non-empty intersection
    m = length(Y.array) # can't use array(Y) ?
    i = 1
    local W
    @inbounds while i <= m
        W = intersection(X, Y.array[i])
        !isempty(W) && break
        i += 1
    end
    if i == m+1
        return EmptySet(dim(Y))
    end

    # add all other non-empty intersections
    out = [W]
    @inbounds for j in i+1:m
        W = intersection(X, Y.array[j])
        !isempty(W) && push!(out, W)
    end
    #println(typeof(overapproximate(UnionSetArray(out), Hyperrectangle)))
    return isempty(out) ? EmptySet(dim(out)) : overapproximate(UnionSetArray(out), Hyperrectangle)
end

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

function _intersection(X::AbstractPolyhedron, Y::AbstractPolyhedron, Z::AbstractPolyhedron, W::AbstractPolyhedron, ::HRepIntersection)
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
# Clustering methods
# =====================================

abstract type AbstractClusteringMethod end

struct NoClustering <: AbstractClusteringMethod
#
end

function cluster(F, idx, ::NoClustering)
    return view(F, idx) # F[idx]
end

struct BoxClustering <: AbstractClusteringMethod
#
end

function cluster(F::Flowpipe{N, ReachSet{N, ST}}, idx, ::BoxClustering) where {N, ST}
    convF = Convexify(view(F, idx))  # Convexify(F[idx])
    return [overapproximate(convF, Hyperrectangle)] # TEMP store as iterable
end

# we use a Zonotope overapproximation of the flowpipe, take thir convex hull, and
# compute its box overapproximation
function cluster(F::Flowpipe{N, TaylorModelReachSet{N}}, idx, method::BoxClustering) where {N, ST}
    Fidx = Flowpipe(view(F, idx))
    charr = [set(overapproximate(Ri, Zonotope)) for Ri in Fidx]
    Xoa = overapproximate(ConvexHullArray(charr), Hyperrectangle)
    return [ReachSet(Xoa, tspan(Fidx))]
end

struct LazyCHClustering <: AbstractClusteringMethod
#
end

function cluster(F, idx, ::LazyCHClustering)
    return [Convexify(view(F, idx))] # Convexify(F[idx])
end

struct ZonotopeClustering <: AbstractClusteringMethod
#
end

# for the generalization to > 2 ses, we iteratively apply the overapprox of the
# CH of two zonotopes, cf LazySets #2154
function cluster(F::Flowpipe{N, RT, VRT}, idx, ::ZonotopeClustering) where {N, ZT<:Zonotope, RT<:ReachSet{N, ZT}, VRT<:AbstractVector{RT}}
    if length(idx) == 1
        return [F[idx]]
    elseif length(idx) == 2
        X = ConvexHull(set(F[idx[1]]), set(F[idx[2]]))
        Y = overapproximate(X, Zonotope) # TODO pass algorithm
        return [ReachSet(Y, tspan(F[idx]))]
    else
        Zaux = overapproximate(ConvexHull(set(F[idx[1]]), set(F[idx[2]])), Zonotope)
        for k in 3:length(idx)
            Zaux = overapproximate(ConvexHull(Zaux, set(F[idx[k]])), Zonotope)
        end
        return [ReachSet(Zaux, tspan(F[idx]))]
    end
end

# convexify and convert to vrep
#C = ReachabilityAnalysis.Convexify(sol[end-aux+1:end])
#Cvertex = convex_hull(vcat([vertices_list(Z) for Z in LazySets.array(set(C))]...)) |> VPolygon

#=
# =====================================
# Methods for checking inclusion
# =====================================

# given two reach-sets A, B check whether f(A) ⊆ g(B) where
# A ⊆ f(A) and B ⊆ g(B)
abstract type AbstractInclusionMethod end

# no pre-processing of the sets
struct LazySetsDefaultInclusion <: AbstractInclusionMethod
#
end
=#
