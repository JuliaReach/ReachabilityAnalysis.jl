# =========================
# Conversion
# =========================

@inline function nodim_msg()
    throw(ArgumentError("to use the `static` option you should pass " *
                        "the system's dimension argument `dim=...`"))
end

# no-op
_reconvert(Ω0::Zonotope{N,Vector{N},Matrix{N}}, static::Val{false}, dim, ngens) where {N} = Ω0
_reconvert(Ω0::Zonotope{N,<:SVector,<:SMatrix}, static::Val{true}, dim) where {N} = Ω0

# convert any zonotope to be represented with regular arrays
function _reconvert(Ω0::Zonotope, static::Val{false}, dim, ngens)
    return Zonotope(Vector(Ω0.center), Matrix(Ω0.generators))
end

# convert any zonotope to be represented with static arrays
function _reconvert(Ω0::Zonotope{N,VN,MN}, static::Val{true}, dim::Missing,
                    ngens::Missing) where {N,VN,MN}
    n, m = size(Ω0.generators)
    return _reconvert(Ω0, static, Val(n), Val(m))
end

function _reconvert(Ω0::Zonotope{N,VN,MN}, static::Val{true}, dim::Val{n},
                    ngens::Val{p}) where {N,VN,MN,n,p}
    G = Ω0.generators
    m = size(G, 2)
    c = SVector{n,N}(Ω0.center)

    if m == p
        return Zonotope(c, SMatrix{n,p}(G))

    elseif m < p
        # extend with zeros
        Gext = hcat(SMatrix{n,m}(G), zeros(MMatrix{n,p - m,N}))
        return Zonotope(c, Gext)

    else
        throw(ArgumentError("can't reconvert a zontope with $m generators to a " *
                            "zonotope with $p generators; you should reduce it first"))
    end
end

# no-op
function _reconvert(Ω0::Hyperrectangle{N,Vector{N},Vector{N}}, static::Val{false},
                    dim::Missing) where {N}
    return Ω0
end
function _reconvert(Ω0::Hyperrectangle{N,Vector{N},Vector{N}}, static::Val{true},
                    dim::Missing) where {N}
    return nodim_msg()
end
_reconvert(Ω0::Hyperrectangle{N,<:SVector,<:SVector}, static::Val{true}, dim) where {N} = Ω0
_reconvert(Ω0::Hyperrectangle{N,<:SVector,<:SVector}, static::Val{true}, dim::Val) where {N} = Ω0

# convert any Hyperrectangle to be represented with regular arrays
function _reconvert(Ω0::Hyperrectangle, static::Val{false}, dim::Missing)
    return Ω0 = Hyperrectangle(Vector(Ω0.center), Matrix(Ω0.radius); check_bounds=false)
end

# convert any Hyperrectangle to be represented with static arrays
function _reconvert(Ω0::Hyperrectangle{N}, static::Val{true}, dim::Val{n}) where {N,n}
    #n = length(Ω0.center) # dimension
    return Ω0 = Hyperrectangle(SVector{n,N}(Ω0.center), SVector{n,N}(Ω0.radius); check_bounds=false)
end

# no-op
_reconvert(Φ::Matrix{N}, static::Val{false}, dim) where {N} = Φ
_reconvert(Φ::IntervalMatrix{N}, static::Val{false}, dim) where {N} = Φ
_reconvert(Φ::AbstractMatrix, static::Val{false}, dim) = Matrix(Φ)
_reconvert(Φ::SMatrix, static::Val{true}, dim) = Φ
_reconvert(Φ::SMatrix, static::Val{true}, dim::Val) = Φ
_reconvert(Φ::SMatrix, static::Val{true}, dim::Missing) = Φ
function _reconvert(Φ::AbstractMatrix, static::Val{true}, dim::Missing)
    return _reconvert(Φ, static, Val(size(Φ, 1)))
end

function _reconvert(Φ::AbstractMatrix{N}, static::Val{true}, dim::Val{n}) where {N,n}
    #n = size(Φ, 1)
    return Φ = SMatrix{n,n,N,n * n}(Φ)
end

function _reconvert(Φ::IntervalMatrix{N,IN,Matrix{IN}}, static::Val{true},
                    dim::Val{n}) where {N,IN,n}
    #n = size(Φ, 1)
    return Φ = IntervalMatrix(SMatrix{n,n,IN,n * n}(Φ))
end

# AbstractVPolygon is not yet available in LazySets
const VPOLY{N,VN} = Union{VPolygon{N,VN},VPolytope{N,VN}}

# no-op
_reconvert(V::VPOLY{N,VN}, static::Val{false}, dim) where {N,VN<:AbstractVector{N}} = V
_reconvert(V::VPOLY{N,VN}, static::Val{true}, dim) where {N,VN<:SVector{N}} = V
_reconvert(V::VPOLY{N,VN}, static::Val{true}, dim::Val) where {N,VN<:SVector{N}} = V
_reconvert(V::VPOLY{N,VN}, static::Val{true}, dim::Missing) where {N,VN<:SVector{N}} = V

function _reconvert(V::VPOLY{N,VN}, static::Val{true},
                    dim::Val{n}) where {N,VN<:AbstractVector{N},n}
    VP = n == 2 ? VPolygon : VPolytope
    return VP([SVector{n,N}(vi) for vi in vertices_list(V)])
end

# dimension is missing
function _reconvert(V::VPOLY{N,VN}, static::Val{true}, dim::Missing) where {N,VN<:AbstractVector{N}}
    return _reconvert(V, static, Val(dim(V)))
end

function Base.convert(::Type{Hyperrectangle{N,Vector{N},Vector{N}}},
                      H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {N,L}
    return Hyperrectangle(Vector(H.center), Vector(H.radius))
end

# duck-typing sampling functions
LazySets._default_sampler(X::IA.Interval) = LazySets._default_sampler(convert(Interval, X))
LazySets._default_sampler(X::IntervalBox) = LazySets._default_sampler(convert(Hyperrectangle, X))
LazySets._default_sampler(X::AbstractVector{<:Real}) = LazySets._default_sampler(Singleton(X))

LazySets.sample(X::IA.Interval, d::Integer; kwargs...) = sample(convert(Interval, X), d; kwargs...)
function LazySets.sample(X::IntervalBox, d::Integer; kwargs...)
    return sample(convert(Hyperrectangle, X), d; kwargs...)
end
function LazySets.sample(X::AbstractVector{<:Real}, d::Integer; kwargs...)
    return sample(Singleton(X), d; kwargs...)
end

LazySets.low(dom::IntervalBox) = inf.(dom.v)
LazySets.low(dom::IntervalBox, i::Int) = inf(dom.v[i])
LazySets.high(dom::IntervalBox) = sup.(dom.v)
LazySets.high(dom::IntervalBox, i::Int) = sup(dom.v[i])

LazySets.low(dom::IA.Interval) = inf(dom)
function LazySets.low(dom::IA.Interval, i::Int)
    @assert i == 1
    return inf(dom)
end
LazySets.high(dom::IA.Interval) = sup(dom)
function LazySets.high(dom::IA.Interval, i::Int)
    @assert i == 1
    return sup(dom)
end

# ------------------------------------------------
# Functions to handle splitting of IntervalBoxes
# TODO refactor to LazySets
# See also: LazySets#2651, IntervalArithmetic#444
# ------------------------------------------------

function Base.convert(HT::Type{Hyperrectangle{N,Vector{N},Vector{N}}},
                      H::AbstractHyperrectangle) where {N}
    c = convert(Vector{N}, LazySets.center(H))
    r = convert(Vector{N}, radius_hyperrectangle(H))
    return Hyperrectangle(c, r)
end

function Base.convert(HT::Type{Hyperrectangle{N,Vector{N},Vector{N}}},
                      B::IntervalBox{D,N}) where {D,N}
    H = convert(Hyperrectangle, B)
    return convert(Hyperrectangle{N,Vector{N},Vector{N}}, H)
end

function Base.split(B::IntervalBox{D,N}, partition::AbstractVector{Int}) where {D,N}
    H = convert(Hyperrectangle{N,Vector{N},Vector{N}}, B)
    return split(H, partition)
end

# =========================
# Projection
# =========================

# extend LazySets concrete projection for other arguments
LazySets.project(X::LazySet, vars::NTuple{D,Int}) where {D} = project(X, collect(vars))
LazySets.project(X::LazySet; vars) = project(X, vars)

# extend LazySets lazy projection for other arguments
LazySets.Projection(X::LazySet, vars::NTuple{D,Int}) where {D} = Projection(X, collect(vars))
LazySets.Projection(X::LazySet; vars) = Projection(X, vars)

# ===============================
# Splitting
# ===============================

# split the symmetric box [-1, 1]^n in nparts in each dimension
function _split_symmetric_box(n::Int, nparts::Int)
    return IA.mince(symBox(n), nparts)
end

# split the symmetric box [-1, 1]^n in partition[i] parts for each dimension
function _split_symmetric_box(D::Int, partition::Vector{Int})
    S = BallInf(zeros(D), 1.0)
    Sp = split(S, partition)
    IB = convert.(IntervalBox, Sp)
    return [IB[i] for i in 1:length(IB)]
end

# ====================================
# Disjointness checks
# ====================================

abstract type AbstractDisjointnessMethod end

# we pass the sets to the disjointness check algorithm without pre-processing
struct FallbackDisjointness <: AbstractDisjointnessMethod end

const NoEnclosure = FallbackDisjointness

# we overapproximate the reach-set with a zonotope, then make the disjointness check
struct ZonotopeEnclosure <: AbstractDisjointnessMethod end

# we overapproximate the reach-set with a hyperrectangle, then make the disjointness check
struct BoxEnclosure <: AbstractDisjointnessMethod end

# this is a dummy disjointness check which returns "false" irrespective of the value of its arguments
struct Dummy <: AbstractDisjointnessMethod end

# --------------------------------------------------------------------
# Methods to evaluate disjointness
# --------------------------------------------------------------------

# fallbacks
_isdisjoint(X::LazySet, Y::LazySet, ::FallbackDisjointness) = isdisjoint(X, Y)
_isdisjoint(X, Y) = _isdisjoint(X, Y, FallbackDisjointness())

# -----------------------------------------------
# Disjointness checks between specific set types
# -----------------------------------------------

# H : {x : ax <= b}, one-dimensional with a != 0
@commutative function _isdisjoint(X::Interval, H::HalfSpace)
    a = H.a[1]
    b = H.b
    N = promote_type(eltype(X), eltype(H))
    if a > zero(N)
        return !_leq(min(X), b / a)
    else
        return !_geq(max(X), b / a)
    end
end

# H : {x : ax = b}, one-dimensional with a != 0
@commutative function _isdisjoint(X::Interval, H::Hyperplane)
    q = H.b / H.a[1]
    return !_geq(q, min(X)) || !_leq(q, max(X))
end

# ====================================================================
# Concrete intersection
#
#  Methods to compute the intersection between two or more sets
# ====================================================================

# -----------------------
# Auxiliary functions
# -----------------------

# converts the normal vector of a list of half-spaces to be a Vector
const VECH{N,VT} = Vector{HalfSpace{N,VT}}
_to_vec(c::HalfSpace{N,Vector{N}}) where {N} = c
_to_vec(c::HalfSpace{N,VT}) where {N,VT<:AbstractVector{N}} = HalfSpace(Vector(c.a), c.b)
_to_vec(x::VECH{N,Vector{N}}) where {N} = x
_to_vec(x::VECH{N,VT}) where {N,VT<:AbstractVector{N}} = broadcast(_to_vec, x)

# concatenates lists of half-spaces such that the normal vectors of the final list
# are all Vector
_vcat(args::VECH...) = vcat(map(_to_vec, args)...)

# ------------------------
# FallbackIntersection
# ------------------------

abstract type AbstractIntersectionMethod end

# "fallback" implementation that uses LazySets intersection(X, Y)
struct FallbackIntersection{T} <: AbstractIntersectionMethod
    backend::T
end

_intersection(X::LazySet, Y::LazySet, ::FallbackIntersection) = intersection(X, Y)
_intersection(X, Y) = _intersection(X, Y, FallbackIntersection())

FallbackIntersection() = FallbackIntersection(nothing)

has_backend(alg::FallbackIntersection) = !isnothing(alg.backend)

function _intersection(X::AbstractPolyhedron{N}, Y::AbstractPolyhedron{N},
                       alg::FallbackIntersection) where {N}
    if has_backend(alg)
        return intersection(X, Y; backend=alg.backend)
    else
        return intersection(X, Y)
    end
end

# ------------------------
# HRepIntersection
# ------------------------

# evaluate X ∩ Y exactly using the constraint representation of X and Y
# evaluate X₁ ∩ ⋯ ∩ Xₖ using the constraint representation of each Xᵢ
#
# TODO Annotate normal vector types?
# struct HRepIntersection{SX, SY} <: AbstractIntersectionMethod end
# setrep(int_method::HRepIntersection{SX<:AbstractPolytope}, SY<:AbstractPolyhedron}) = SX
#
struct HRepIntersection <: AbstractIntersectionMethod
    #
end

setrep(::HRepIntersection) = HPolytope{Float64,Vector{Float64}}

function _intersection(X::AbstractPolyhedron, Y::AbstractPolyhedron, ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    out = _vcat(clist_X, clist_Y)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

function _intersection(X::AbstractPolyhedron, Y::AbstractPolyhedron, Z::AbstractPolyhedron,
                       ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    clist_Z = constraints_list(Z)
    out = _vcat(clist_X, clist_Y, clist_Z)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

function _intersection(X::LazySet, Y::AbstractPolyhedron, Z::AbstractPolyhedron,
                       W::AbstractPolyhedron, ::HRepIntersection)
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)
    clist_Z = constraints_list(Z)
    clist_W = constraints_list(W)
    out = _vcat(clist_X, clist_Y, clist_Z, clist_W)
    success = remove_redundant_constraints!(out)
    return (success, out)
end

# ------------------------
# BoxIntersection
# ------------------------

struct BoxIntersection <: AbstractIntersectionMethod
    #
end

setrep(::BoxIntersection) = Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}

# ------------------------
# TemplateHullIntersection
# ------------------------

# evaluate X ∩ Y approximately using support function evaluations
#
# if lazy = false (default) use the property that the support function
# of X ∩ Y along direction d is not greater
# ρ(d, X ∩ Y) <= min(ρ(d, X), ρ(d, Y))
# by the same token, compute X₁ ∩ ⋯ ∩ Xₖ approximately using the same property
# if the template is provided, we have TN<:AbstractDirections{N, VN}
# otherwise, the constraints of X and Y are used and TN is Missing
#
# if lazy = true, use specialized approximation of lazy intersections
# assuming that Y is polyhedral
struct TemplateHullIntersection{N,VN,TN,L} <: AbstractIntersectionMethod
    dirs::TN
    lazy::L
end

# constructor with template directions provided
function TemplateHullIntersection(dirs::TN; lazy=false) where {N,VN,TN<:AbstractDirections{N,VN}}
    lazy_val = Val(lazy)
    return TemplateHullIntersection{N,VN,TN,typeof(lazy_val)}(dirs, lazy_val)
end

# constructor without template directions => directions are missing until evaluated
function TemplateHullIntersection{N,VN}(; lazy=false) where {N,VN<:AbstractVector{N}}
    lazy_val = Val(lazy)
    return TemplateHullIntersection{N,VN,Missing,typeof(lazy_val)}(missing, lazy_val)
end

function TemplateHullIntersection(; lazy=false)
    return TemplateHullIntersection{Float64,Vector{Float64}}(; lazy=lazy)
end

setrep(::TemplateHullIntersection{N,VN}) where {N,VN} = HPolytope{N,VN}
function setrep(::TemplateHullIntersection{N,SEV}) where {N,SEV<:SingleEntryVector{N}}
    return Union{HPolytope{N,SEV},HPolytope{N,Vector{N}}}
end
function setrep(::TemplateHullIntersection{N,SP}) where {N,SP<:SparseVector{N}}
    return Union{HPolytope{N,SP},HPolytope{N,Vector{N}}}
end

# if the template directions is missing => use the constraints of both X and Y
# doesn't remove redundant constraints
function _intersection(X::LazySet, Y::LazySet,
                       method::TemplateHullIntersection{N,VN,Missing,Val{false}}) where {N,VN}
    clist_X = constraints_list(X)
    clist_Y = constraints_list(Y)

    out = Vector{HalfSpace{N,VN}}()
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
function _intersection(X::LazySet, Y::LazySet,
                       method::TemplateHullIntersection{N,VN,TN,Val{false}}) where {N,VN,
                                                                                    TN<:AbstractDirections{N,
                                                                                                           VN}}
    dirs = method.dirs
    out = Vector{HalfSpace{N,VN}}(undef, length(dirs))
    @inbounds for (i, d) in enumerate(dirs)
        d = convert(VN, d)
        b = min(ρ(d, X), ρ(d, Y))
        out[i] = HalfSpace(d, b)
    end
    return HPolytope(out)
end

# compute the min of X ∩ Hi for each Hi in Y (assuming the second set is polyhedral,
# requires the list of constraints of Y) for each template direction d
function _intersection(X::LazySet, Y::LazySet,
                       method::TemplateHullIntersection{N,VN,TN,Val{true}}) where {N,VN,
                                                                                   TN<:AbstractDirections{N,
                                                                                                          VN}}
    dirs = method.dirs
    out = Vector{HalfSpace{N,VN}}(undef, length(dirs))
    H = constraints_list(Y)

    @inbounds for (i, d) in enumerate(dirs)
        d = convert(VN, d)
        b = minimum(ρ(d, X ∩ H[j]) for j in eachindex(H))
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
struct FallbackInclusion <: AbstractInclusionMethod
    #
end

_iscontained(X::LazySet, Y::LazySet, ::FallbackInclusion) = LazySets .⊆ (X, Y)
_iscontained(X, Y) = _iscontained(X, Y, FallbackInclusion())
