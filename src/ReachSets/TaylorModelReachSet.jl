# ================================================================
# Taylor model reach set
# ================================================================

using TaylorModels: TaylorModel1, TaylorN, fp_rpa

"""
    TaylorModelReachSet{N, S} <: AbstractTaylorModelReachSet{N}

Reach-set representation consisting of a vector of taylor models in one variable
(the "time" variable) whose coefficients are multivariate polynomials
(the "space" variables).

### Notes

The parameter `N` refers to the numerical type of the representation.
The space variables are assumed to be normalized to the interval `[-1, 1]`.
It is assumed that the time domain is the same for all components.
"""
struct TaylorModelReachSet{N, S} <: AbstractTaylorModelReachSet{N}
    X::Vector{TaylorModel1{TaylorN{S}, N}}
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

# constructor with a time point
function TaylorModelReachSet(X::Vector{TaylorModel1{TaylorN{S}, N}}, t::Real) where {S, N}
    return TaylorModelReachSet(X, interval(t))
end

# ======================
# Intersection methods
# ======================

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

# ======================
# Inclusion checks
# ======================

function Base.:⊆(R::TaylorModelReachSet, X::LazySet)
    Rz = overapproximate(R, Zonotope)
    return set(Rz) ⊆ X
end

# ======================
# Disjointness checks
# ======================

_is_intersection_empty(R::TaylorModelReachSet, Y::LazySet, method::FallbackDisjointness) = isdisjoint(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::LazySet, method::ZonotopeEnclosure) = isdisjoint(set(overapproximate(R, Zonotope)), Y)

# FIXME when UnionSet, UnionSetArray <: LazySet
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::FallbackDisjointness) = isdisjoint(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::ZonotopeEnclosure) = isdisjoint(set(overapproximate(R, Zonotope)), Y)
_is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet, UnionSetArray}, method::BoxEnclosure) = isdisjoint(set(overapproximate(R, Hyperrectangle)), Y)

_is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet, method::FallbackDisjointness) = isdisjoint(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))
_is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet, method::ZonotopeEnclosure) = isdisjoint(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))


@commutative function _is_intersection_empty(R::TaylorModelReachSet, Y::Union{LazySet, UnionSet, UnionSetArray}, method::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return isdisjoint(set(Z), Y)
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

# ==============================================================
# Conversion and overapproximation of Taylor model reach-sets
# ==============================================================

# no-op
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:TaylorModelReachSet}) where {N}
    return R
end

function overapproximate(R::TaylorModelReachSet{N}, ::BoxEnclosure) where {N}
    return overapproximate(R, Hyperrectangle)
end

function overapproximate(R::TaylorModelReachSet{N}, ::ZonotopeEnclosure) where {N}
    return overapproximate(R, Zonotope)
end

# tdom defined in the methods below is the same as Δt - t0, but the domain inclusion check
# in TM.evauate may fail, so we unpack the domain again here; also note that
# by construction the TMs in time are centered at zero

# overapproximate taylor model reachset with one hyperrectangle
function _overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle}) where {N}
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

LazySets.box_approximation(R::TaylorModelReachSet) = _overapproximate(R, Hyperrectangle)

# overapproximate taylor model reachset with several hyperrectangles,
# by splitting in time (ntdiv), and/or space (with either nsdiv for uniform partitions
# or by specifying a partition (vector of integers))
# e.g. if the partition is uniform for each dimension, the number of returned sets
# is nsdiv^D * ntdiv
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Hyperrectangle};
                         partition=nothing, nsdiv=nothing, ntdiv=1) where {N}

    if !isnothing(partition) && !isnothing(nsdiv)
        throw(ArgumentError("either `partition` or `nsdiv` should be specified, not both"))
    end
    if isnothing(partition) && isnothing(nsdiv)
        nsdiv = 1
    end

    # no splitting
    if isnothing(partition) && nsdiv == 1 && ntdiv == 1
        return _overapproximate(R, Hyperrectangle)
    end

    D = dim(R)
    X = set(R)

    # time domain splittig
    tspdiv = IntervalArithmetic.mince(tspan(R), ntdiv)
    tdom = domain(R)
    Δtdiv = IntervalArithmetic.mince(tdom, ntdiv)

    # evaluate the Taylor model in (normalized) time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = [evaluate(X, Δti) for Δti in Δtdiv]

    # partition the symmetric box
    S = symBox(D)
    if isnothing(partition)
        if nsdiv == 1
            partition = [S]
        else
            # FIXME (also below) may use IntervalArithmetic.mince directly
            partition = fill(nsdiv, D)
            Sdiv = LazySets.split(convert(Hyperrectangle, S), partition)
            partition = convert.(IntervalBox, Sdiv)
        end

    else
        Sdiv = LazySets.split(convert(Hyperrectangle, S), partition)
        partition = convert.(IntervalBox, Sdiv)
    end
    nparts = length(partition)

    # preallocate output vectors
    ST = Hyperrectangle{N, SVector{D, N}, SVector{D, N}}
    X̂out = Vector{ST}(undef, nparts * ntdiv)
    R̂out = Vector{ReachSet{N, ST}}(undef, nparts * ntdiv)

    # evaluate the spatial variables in the symmetric box
    @inbounds for k in 1:ntdiv
        for j in 1:nparts
            Bj = partition[j]
            Xk = X_Δt[k]
            X̂ib = IntervalBox([evaluate(Xk[i], Bj) for i in 1:D])

            idx = j + (k-1)*nparts
            X̂out[idx] = convert(Hyperrectangle, X̂ib)

            R̂out[idx] = ReachSet(X̂out[idx], tspdiv[k])
        end
    end

    return R̂out
end

# overapproximate taylor model reachset with one zonotope
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope};
                         Δt::TimeInterval=tspan(R),
                         dom::IntervalBox=symBox(dim(R))) where {N}
    @assert dim(R) == dim(dom) "the dimension of the reach-set should match the dimension of the domain, but they are $(dim(R)) and $(dim(dom)) respectively"
    if !(dom ⊆ symBox(dim(R)))
        throw(ArgumentError("`dom` must be a subset of [-1, 1]^n"))
    end
    if !(Δt ⊆ tspan(R))
        throw(ArgumentError("the given time range $Δt does not belong to the " *
                            "reach-set's time span, $(tspan(R))"))
    end

    # dimension of the reachset
    n = dim(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X = set(R)
    tdom = Δt - tstart(R)  # normalize time (to TM-internal time)
    tdom = tdom ∩ domain(R)  # intersection handles round-off errors
    X_Δt = evaluate(X, tdom)

    # builds the associated taylor model for each coordinate j = 1...n
    #  X̂ is a TaylorModelN whose coefficients are intervals
    X̂ = [TaylorModelN(X_Δt[j], zeroI, zeroBox(n), dom) for j in 1:n]

    # compute floating point rigorous polynomial approximation
    # fX̂ is a TaylorModelN whose coefficients are floats
    fX̂ = fp_rpa.(X̂)

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
        X̂ib = [TaylorModelN(X_Δt[j], zeroI, x0, Bi) for j in 1:D]
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
    tn = t - tstart(R)  # normalize time (to TM-internal time)
    tn = interval(tn) ∩ domain(R)  # intersection handles round-off errors
    X_Δt = evaluate(X, tn)
    n = dim(R)
    X̂ = [TaylorModelN(X_Δt[j], zeroI, zeroBox(n), symBox(n)) for j in 1:n]
    fX̂ = fp_rpa.(X̂)
    Zi = overapproximate(fX̂, Zonotope, remove_zero_generators=remove_zero_generators)
    Δt = TimeInterval(t, t)
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

# ======================
# Remainder handling
# ======================

using TaylorModels: shrink_wrapping!

function _shrink_wrapping(R::TaylorModelReachSet)
    rem = remainder(R)
    dt = domain(R)
    n = dim(R)

    # if all remainders are zero => nothing to do
    all(iszero, rem) && return R

    # transform into a TaylorN whose coeffs are floats
    X = set(R)
    Xev = evaluate.(X, dt)
    W = [TaylorModelN(Xev[i], zeroI, zeroBox(n), symBox(n)) for i in 1:n]
    Wfp = fp_rpa.(W)

    # absorb remainder in the polynomial part
    shrink_wrapping!(Wfp)

    # transform back to a TaylorModel1 of TaylorN
    orderT = get_order(set(R)[1])
    p = [Taylor1(TaylorN(polynomial(Wfp[i])), orderT) for i in 1:n]
    Y = [TaylorModel1(p[i], zeroI, zeroI, dt) for i in 1:n]

    return TaylorModelReachSet(Y, tspan(R))
end
