# ================================================================
# Taylor model reach set
# ================================================================

using TaylorModels: TaylorModel1, TaylorN, fp_rpa

"""
    TaylorModelReachSet{N} <: AbstractTaylorModelReachSet{N}

Taylor model reach-set represented as a vector of taylor models in one variable
(namely, the "time" variable) whose coefficients are multivariate polynomials
(namely in the "space" variables). It is assumed that the time domain is the same
for all components.

### Notes

The parameter `N` refers to the numerical type of the representation.

In `TMJets`, the space variables are normalized to the interval `[-1, 1]`.
"""
struct TaylorModelReachSet{N} <: AbstractTaylorModelReachSet{N}
    X::Vector{TaylorModel1{TaylorN{N}, N}}
    #t0::Float64
    #δt::Float64
    Δt::TimeInterval
end

# interface functions
@inline set(R::TaylorModelReachSet) = R.X
@inline setrep(::TaylorModelReachSet{N}) where {N} = Vector{TaylorModel1{TaylorN{N}, N}}
@inline setrep(::Type{TaylorModelReachSet{N}}) where {N} = Vector{TaylorModel1{TaylorN{N}, N}}
@inline tstart(R::TaylorModelReachSet) = tstart(R.Δt) # t0 - δt + Δt0.lo
@inline tend(R::TaylorModelReachSet) = tend(R.Δt) # t0 + Δt0.hi
@inline tspan(R::TaylorModelReachSet) = R.Δt # # Interval(tstart(R), tend(R))
@inline dim(R::TaylorModelReachSet) = get_numvars()
@inline vars(R::TaylorModelReachSet) = Tuple(Base.OneTo(length(R.X)),)

# overload getter functions for the taylor model
# we assume that the first element is representative
domain(R::TaylorModelReachSet) = domain(first(R.X)) # normalized time domain
remainder(R::TaylorModelReachSet) = remainder.(R.X)
polynomial(R::TaylorModelReachSet) = polynomial.(R.X)
get_order(R::TaylorModelReachSet) = get_order.(R.X)
expansion_point(R::TaylorModelReachSet) = [Xi.x0 for Xi in R.X]

function shift(R::TaylorModelReachSet, t0::Number)
    return TaylorModelReachSet(set(R), tspan(R) + t0)
end

function reconstruct(R::TaylorModelReachSet, X)
    return TaylorModelReachSet(X, tspan(R))
end

# TODO remove but add comment in docs
function project(R::TaylorModelReachSet, vars::NTuple{D, M}) where {D, M<:Integer}
    throw(ArgumentError("the concrete projection of Taylor model reach-set is not " *
            "available; try first to overapproximate the Taylor model and the  project"))
end

function overapproximate(R::TaylorModelReachSet{N}, ::BoxEnclosure) where {N}
    return overapproximate(R, Hyperrectangle)
end

function overapproximate(R::TaylorModelReachSet{N}, ::ZonotopeEnclosure) where {N}
    return overapproximate(R, Zonotope)
end

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
    return isempty(out) ? EmptySet(dim(out)) : overapproximate(UnionSetArray(out), Hyperrectangle)
end

# ==============================================================
# Conversion and overapproximation of Taylor model reach-sets
# ==============================================================

# no-op
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:TaylorModelReachSet}) where {N}
    return R
end

# tdom defined in the methods below is the same as Δt - t0, but the domain inclusion check
# in TM.evauate may fail, so we unpack the domain again here; also note that
# by construction the TMs in time are centered at zero

# overapproximate taylor model reachset with one hyperrectangle
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle}) where {N}
    # dimension of the reachset
    D = dim(R)

    # normalized time domain
    tdom = domain(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = evaluate(set(R), tdom)

    # evaluate the spatial variables in the symmetric box
    Bn = symBox(D)
    X̂ib = IntervalBox([evaluate(X_Δt[i], Bn) for i in 1:D]...)
    X̂ = convert(Hyperrectangle, X̂ib)

    Δt = tspan(R)
    return ReachSet(X̂, Δt)
end

LazySets.box_approximation(R::TaylorModelReachSet) = overapproximate(R, Hyperrectangle)

# overapproximate taylor model reachset with several hyperrectangles
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle}, nparts) where {N}
    # dimension of the reachset
    D = dim(R)

    # normalized time domain
    tdom = domain(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = evaluate(set(R), tdom)

    # evaluate the spatial variables in the symmetric box
    partition = IA.mince(symBox(D), nparts)
    X̂ = Vector{Hyperrectangle{N, SVector{D, N}, SVector{D, N}}}(undef, length(partition))
    @inbounds for (i, Bi) in enumerate(partition)
        X̂ib = IntervalBox([evaluate(X_Δt[i], Bi) for i in 1:D])
        X̂[i] = convert(Hyperrectangle, X̂ib)
    end

    Δt = tspan(R)
    #return ReachSet(UnionSetArray(X̂), Δt) # but UnionSetArray is not yet a lazyset
    return ReachSet(ConvexHullArray(X̂), Δt)
end

# overapproximate taylor model reachset with one zonotope
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope}) where {N}
    # dimension of the reachset
    n = dim(R)

    # normalized time domain
    tdom = domain(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X = set(R)
    X_Δt = evaluate(X, tdom)

    # builds the associated taylor model for each coordinate j = 1...n
    #  X̂ is a TaylorModelN whose coefficients are intervals
    X̂ = [TaylorModelN(X_Δt[j], zeroI, zeroBox(n), symBox(n)) for j in 1:n]

    # compute floating point rigorous polynomial approximation
    # fX̂ is a TaylorModelN whose coefficients are floats
    fX̂ = fp_rpa.(X̂)

    Δt = tspan(R)
    # LazySets can overapproximate a Taylor model with a Zonotope
    Zi = overapproximate(fX̂, Zonotope)
    return ReachSet(Zi, Δt)
end

# overapproximate taylor model reachset with several zonotopes given the partition,
# that specifies the way in which the symmetric box [-1, 1]^n is split:
#
# - if `partition` is an integer, makes as uniform partition for each coordinate
# - if `partition` is a vector of integers, split the domain according to each element in partition
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope}, partition::Union{Int, Vector{Int}}) where {N}
    # dimension of the reachsets
    D = dim(R)

    # normalized time domain
    tdom = domain(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X = set(R)
    X_Δt = evaluate(X, tdom)

    # evaluate the spatial variables in the symmetric box
    part = _split_symmetric_box(D, partition)
    fX̂ = Vector{Vector{TaylorModelN{length(X_Δt), N, N}}}(undef, length(part))
    @inbounds for (i, Bi) in enumerate(part)
        x0 = IntervalBox(mid.(Bi))
        X̂ib = [TaylorModelN(X_Δt[j], X[j].rem, x0, Bi) for j in 1:D]
        fX̂[i] = fp_rpa.(X̂ib)
    end
    Z = overapproximate.(fX̂, Zonotope)
    Δt = tspan(R)
    #return ReachSet(UnionSetArray(Z), Δt) # but UnionSetArray is not yet a lazyset
    return ReachSet(ConvexHullArray(Z), Δt)
end

# evaluate at a given time and overapproximate the resulting set with a zonotope
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope}, t::AbstractFloat;
                         remove_zero_generators=true) where {N}
    @assert t ∈ tspan(R) "the given time point $t does not belong to the reach-set's time span, $(tspan(R))"

    X = set(R)
    Δtn = domain(R)
    if t == tstart(R)
        tn = zero(N)
    elseif t == tend(R)
        tn = Δtn.hi
    else
        tn = t - tstart(R)
    end
    X_Δt = evaluate(X, tn)
    n = dim(R)
    X̂ = [TaylorModelN(X_Δt[j], X[j].rem, zeroBox(n), symBox(n)) for j in 1:n]
    fX̂ = fp_rpa.(X̂)
    Zi = overapproximate(fX̂, Zonotope, remove_zero_generators=remove_zero_generators)
    Δt = TimeInterval(t, t)
    return ReachSet(Zi, Δt)
end

function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope}, Δt::TimeInterval) where {N}
    @assert Δt ⊆ tspan(R) "the given time point $t does not belong to the reach-set's time span, $(tspan(R))"

    X = set(R)
    Δtn = domain(R)

    if Δt.lo == tstart(R)
        dtn_lo = zero(N)
    elseif Δt.lo == tend(R)
        dtn_lo = Δtn.hi
    else
        dtn_lo = Δt.lo - tstart(R)
    end

    if Δt.hi == tstart(R)
        dtn_hi = zero(N)
    elseif Δt.hi == tend(R)
        dtn_hi = Δtn.hi
    else
        dtn_hi = Δt - tstart(R)
    end

    dtn = IA.Interval(dtn_lo, dtn_hi)
    X_Δt = evaluate(X, dtn)
    n = dim(R)
    X̂ = [TaylorModelN(X_Δt[j], X[j].rem, zeroBox(n), symBox(n)) for j in 1:n]
    fX̂ = fp_rpa.(X̂)
    Zi = overapproximate(fX̂, Zonotope)

    return ReachSet(Zi, Δt)
end

# convert a hyperrectangular set to a taylor model reachset
function convert(::Type{<:TaylorModelReachSet}, H::AbstractHyperrectangle{N};
                 orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}

    n = dim(H)
    x = set_variables("x", numvars=n, order=orderQ)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N}, N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # for each variable i = 1, .., n, compute the linear polynomial that covers
    # the line segment corresponding to the i-th edge of H
    @inbounds for i in 1:n
        α = radius_hyperrectangle(H, i)
        β = center(H, i)
        pi = α * x[i] + β
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), zeroI, zeroI, Δtn)
    end

    return TaylorModelReachSet(vTM, Δt)
end

# overapproximate a hyperrectangular set with a taylor model reachset, fallback to convert
function overapproximate(H::AbstractHyperrectangle{N}, T::Type{<:TaylorModelReachSet};
                         orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}
    convert(T, H; orderQ=orderQ, orderT=orderT, Δt=Δt)
end

# overapproximate a zonotopic set with a taylor model reachset
# FIXME pass algorithm option to choose between using parallelotope oa or  order reduction
function overapproximate(Z::AbstractZonotope{N}, ::Type{<:TaylorModelReachSet};
                         orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI,
                         indices=1:dim(Z), box_reduction=false) where {N}

    n = dim(Z)
    x = set_variables("x", numvars=n, order=orderQ)

    if order(Z) > 1
        if box_reduction
            # diagonal generators matrix
            Z = _reduce_order(Z, 1)
        else
            # indices selects the indices that we want to keep
            Z = _overapproximate_hparallelotope(Z, indices)
        end
    end
    c = center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N}, N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # for each variable i = 1, .., n, compute the linear polynomial that covers
    # the line segment corresponding to the i-th edge of Z
    @inbounds for i in 1:n
        pi = c[i] + sum(view(G, i, :) .* x)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), zeroI, zeroI, Δtn)
    end

    return TaylorModelReachSet(vTM, Δt)
end

# convert a zonotopic set to a taylor model reachset (only if it has order 1)
function convert(TM::Type{<:TaylorModelReachSet}, Z::AbstractZonotope; kwargs...)
    if order(Z) == 1
        return overapproximate(Z, TM; kwargs...)

    else
        throw(ArgumentError("can't convert a zonotope of order $(order(Z)) to a `TaylorModelReachSet`; " *
                            "use `overapproximate(Z, TaylorModelReachSet)` instead"))
    end
end

function convert(TM::Type{<:TaylorModelReachSet}, R::AbstractLazyReachSet; kwargs...)
    return convert(TM, set(R); Δt=tspan(R), kwargs...)
end

function overapproximate(R::AbstractLazyReachSet, TM::Type{<:TaylorModelReachSet}; kwargs...)
    return overapproximate(set(R), TM; Δt=tspan(R), kwargs...)
end

# zonotope of order 2 with generators matrix G = [M; D] where M is n x n and D is n x n and diagonal
function _overapproximate_structured(Z::AbstractZonotope{N}, ::Type{<:TaylorModelReachSet};
                                     orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}
    n = dim(Z)
    x = set_variables("x", numvars=n, order=orderQ)

    # check structure
    order(Z) == 2 || throw(ArgumentError("this function requires that the order of the zonotope is 2, got $(order(Z))"))

    c = LazySets.center(Z)
    G = genmat(Z)

    M = view(G, :, 1:n)
    D = view(G, :, n+1:2n)
    isdiag(D) || throw(ArgumentError("the columns $(n+1) to $(2n) of the generators matrix do not form a diagonal matrix"))

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N}, N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # for each variable i = 1, .., n, compute the linear polynomial that covers
    # the line segment corresponding to the i-th edge of Z
    @inbounds for i in 1:n
        pi = c[i] + sum(view(M, i, :) .* x)
        di = abs(D[i, i])
        rem = interval(-di, di)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end

    return TaylorModelReachSet(vTM, Δt)
end

function _overapproximate_structured(Zcp::CartesianProduct{N, <:Zonotope, <:Interval}, ::Type{<:TaylorModelReachSet};
                                     orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}
    n = dim(Zcp)
    x = set_variables("x", numvars=n, order=orderQ)

    # check structure
    Z = Zcp.X
    Y = Zcp.Y
    if order(Z) > 1
        Z = remove_redundant_generators(Z)
    end
    order(Z) == 1 || throw(ArgumentError("this function requires that the order of the (lower) zonotope is 1, got $(order(Z))"))

    c = LazySets.center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N}, N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    xstate = x[1:n-1]
    @inbounds for i in 1:n-1
        pi = c[i] + sum(view(G, i, :) .* xstate) + zero(TaylorN(n, order=orderQ))
        rem = interval(0)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end
    # diagonal part
    @inbounds begin
        pi = mid(Y.dat) + zero(TaylorN(1, order=orderQ))
        d = diam(Y.dat) / 2
        rem = interval(-d, d)
        vTM[n] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end
    return TaylorModelReachSet(vTM, Δt)
end

# we assume that if n = size(G, 1), then:
# (size(G) == (n, 2n + 1)) && isdiag(view(G, :, (n+2):(2n+1)))
function _overapproximate_structured_full(Zcp::CartesianProduct{N, <:Zonotope, <:Interval}, ::Type{<:TaylorModelReachSet};
                                          orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}
    n = dim(Zcp) - 1
    x = set_variables("x", numvars=n+1, order=orderQ)

    # check structure
    # not checking structure
    Z = Zcp.X
    c = LazySets.center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N}, N}}(undef, n+1)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # fill rows corresponding to the "zonotope" variables: 1 to nth-variables
    @inbounds for i in 1:n
        pi = c[i] + sum(view(G, i, 1:(n+1)) .* x) + zero(TaylorN(n+1, order=orderQ))
        d = abs(G[i, n + 1 + i])
        rem = interval(-d, d)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end

    # fill the final row, which correspponds to the "interval" variable: (n+1)-th
    I = Zcp.Y.dat
    pi = mid(I) + zero(TaylorN(n+1, order=orderQ))
    d = diam(I) / 2
    rem = interval(-d, d)
    @inbounds vTM[n+1] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)

    return TaylorModelReachSet(vTM, Δt)
end

_is_intersection_empty(R::TaylorModelReachSet, Y::LazySet, method::FallbackDisjointness) = is_intersection_empty(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::LazySet, method::ZonotopeEnclosure) = is_intersection_empty(set(overapproximate(R, Zonotope)), Y)

# FIXME when UnionSet, UnionSetArray <: LazySet
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::FallbackDisjointness) = is_intersection_empty(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::ZonotopeEnclosure) = is_intersection_empty(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::BoxEnclosure) = is_intersection_empty(set(overapproximate(R, Hyperrectangle)), Y)

_is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet, method::FallbackDisjointness) = is_intersection_empty(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))
_is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet, method::ZonotopeEnclosure) = is_intersection_empty(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))


@commutative function _is_intersection_empty(R::TaylorModelReachSet, Y::Union{LazySet, UnionSet, UnionSetArray}, method::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return is_intersection_empty(set(Z), Y)
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe, method::AbstractDisjointnessMethod)
    return false
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe, method::ZonotopeEnclosure)
    return false
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe, method::BoxEnclosure)
    return false
end
