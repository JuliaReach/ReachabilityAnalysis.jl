# ================================================================
# Taylor model reach set
# ================================================================

import ..Overapproximate: _overapproximate, _overapproximate_hparallelotope

"""
    TaylorModelReachSet{N, S} <: AbstractTaylorModelReachSet{N}

Reach-set representation consisting of a vector of taylor models in one variable
(the "time" variable) whose coefficients are multivariate polynomials
(the "space" variables).

### Notes

The parameters `N` and `S` refer to the numerical type of the representation (used for the Taylor model
in time and Taylor series polynomials respectively).
The space variables are assumed to be normalized to the interval `[-1, 1]`.
It is assumed that the time domain is the same for all components.
"""
struct TaylorModelReachSet{N,S} <: AbstractTaylorModelReachSet{N}
    X::Vector{TaylorModel1{TaylorN{S},N}}
    Δt::TimeInterval
end

# interface functions
set(R::TaylorModelReachSet) = R.X
setrep(::TaylorModelReachSet{N}) where {N} = Vector{TaylorModel1{TaylorN{N},N}}
setrep(::Type{TaylorModelReachSet{N,N}}) where {N} = TaylorModelReachSet{N,N}
function setrep(::Type{TaylorModelReachSet{N,IA.Interval{N}}}) where {N}
    return TaylorModelReachSet{N,IA.Interval{N}}
end
setrep(::Type{TaylorModelReachSet{N}}) where {N} = Vector{TaylorModel1{TaylorN{N},N}}
tstart(R::TaylorModelReachSet) = tstart(R.Δt) # t0 - δt + Δt0.lo
tend(R::TaylorModelReachSet) = tend(R.Δt) # t0 + Δt0.hi
tspan(R::TaylorModelReachSet) = R.Δt # # Interval(tstart(R), tend(R))
dim(R::TaylorModelReachSet) = get_numvars()
vars(R::TaylorModelReachSet) = Tuple(Base.OneTo(length(R.X)))

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
function project(::TaylorModelReachSet, ::NTuple{D,Int}) where {D}
    throw(ArgumentError("the concrete projection of a Taylor model reach-set is not " *
                        "available; try to first overapproximate the Taylor model and then project"))
end

# constructor with a time point
function TaylorModelReachSet(X::Vector{TaylorModel1{TaylorN{S},N}}, t::Real) where {S,N}
    return TaylorModelReachSet(X, interval(t))
end

# ======================
# Intersection methods
# ======================

# for TMs, we overapproximate with a zonotope
function _intersection(R::TaylorModelReachSet, Y::LazySet)
    return _intersection(R, Y, FallbackIntersection())
end

function _intersection(R::TaylorModelReachSet, Y::LazySet, ::FallbackIntersection)
    X = set(overapproximate(R, Zonotope))
    return intersection(X, Y)
end

# TODO overapprox Y with a box if it's bounded?
function _intersection(R::TaylorModelReachSet, Y::LazySet, ::BoxIntersection)
    X = set(overapproximate(R, Hyperrectangle))
    out = intersection(X, Y)
    return isempty(out) ? EmptySet(dim(out)) : overapproximate(out, Hyperrectangle)
end

# TODO refactor? See LazySets#2158
function _intersection(R::TaylorModelReachSet, Y::UnionSetArray{N,HT},
                       ::BoxIntersection) where {N,HT<:HalfSpace{N}}
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
    if i == m + 1
        return EmptySet(dim(Y))
    end

    # add all other non-empty intersections
    out = [W]
    @inbounds for j in (i + 1):m
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

function _is_intersection_empty(R::TaylorModelReachSet, Y::LazySet, method::FallbackDisjointness)
    return isdisjoint(set(overapproximate(R, Zonotope)), Y)
end

# TODO when UnionSet, UnionSetArray <: LazySet
function _is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet,UnionSetArray},
                                method::FallbackDisjointness)
    return isdisjoint(set(overapproximate(R, Zonotope)), Y)
end
function _is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet,UnionSetArray},
                                method::ZonotopeEnclosure)
    return isdisjoint(set(overapproximate(R, Zonotope)), Y)
end
function _is_intersection_empty(R::TaylorModelReachSet, Y::Union{UnionSet,UnionSetArray},
                                method::BoxEnclosure)
    return isdisjoint(set(overapproximate(R, Hyperrectangle)), Y)
end

function _is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet,
                                method::FallbackDisjointness)
    return isdisjoint(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))
end
function _is_intersection_empty(R1::TaylorModelReachSet, R2::TaylorModelReachSet,
                                method::ZonotopeEnclosure)
    return isdisjoint(set(overapproximate(R1, Zonotope)), set(overapproximate(R2, Zonotope)))
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet,
                                             Y::Union{LazySet,UnionSet,UnionSetArray},
                                             method::ZonotopeEnclosure)
    Z = overapproximate(R, Zonotope)
    return isdisjoint(set(Z), Y)
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe,
                                             method::AbstractDisjointnessMethod)
    return false
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe,
                                             method::ZonotopeEnclosure)
    return false
end

@commutative function _is_intersection_empty(R::TaylorModelReachSet, X::Universe,
                                             method::BoxEnclosure)
    return false
end

# ======================
# Remainder handling
# ======================

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

# ======================
# Domain handling
# ======================

# Given a vector of multivariate polynomials in the normalized symmetric box
# `[-1, 1]^n`, and a (sub)-domain `dom ⊆ D`, with low/high bounds
# `a, b ∈ R^n` respectively, apply the transformation
# `x[i] <- (a[i] + b[i]) / 2 + (b[i] - a[i]) / 2 * x[i]` for each `i = 1, .., n`,
# which amounts to rescaling and shifting the polynomials according to `dom`.
function _taylor_shift(X::Vector{TaylorN{S}}, dom) where {S}
    x = get_variables()
    n = length(x) # number of variables
    @assert n == length(X) == dim(dom)
    (dom == symBox(n)) && return X

    a, b = low(dom), high(dom)
    transf = [(a[i] + b[i]) / 2 + (b[i] - a[i]) / 2 * x[i] for i in 1:n]
    return [evaluate(X[i], transf) for i in 1:n]
end

# =================================
# Evaluation
# =================================

function evaluate(R::TaylorModelReachSet, Δt::TimeInterval)
    n = dim(R)
    X = set(R)
    Δtn = (Δt - tstart(R)) ∩ domain(R)
    return [fp_rpa(TaylorModelN(evaluate(X[i], Δtn), zeroI, zeroBox(n), symBox(n))) for i in 1:n]
end

evaluate(R::TaylorModelReachSet, t::Real) = evaluate(R, interval(t))

# =================================
# Conversion and overapproximation
# =================================

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
# in TM.evaluate may fail, so we unpack the domain again here; also note that
# by construction the TMs in time are centered at zero

# overapproximate Taylor model reachset with a set
function _overapproximate(R::TaylorModelReachSet, T::Type{<:LazySet};
                          Δt::TimeInterval=tspan(R), dom=symBox(dim(R)), kwargs...)
    # dimension of the reachset
    n = dim(R)

    # consistency checks
    n == dim(dom) ||
        throw(ArgumentError("the dimension of the reach-set should match the dimension of the domain, " *
                            "but they are $(dim(R)) and $(dim(dom)) respectively"))
    dom ⊆ symBox(dim(R)) || throw(ArgumentError("`dom` must be a subset of [-1, 1]^n"))
    Δt ⊆ tspan(R) || throw(ArgumentError("the given time range $Δt does not belong to the " *
                                         "reach-set's time span, $(tspan(R))"))

    # evaluate the Taylor model in time
    X = set(R)
    tdom = Δt - tstart(R)  # normalize time (to TM-internal time)
    tdom = tdom ∩ domain(R)  # intersection handles round-off errors
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = evaluate(X, tdom)

    # transform the domain if needed
    X_Δt = _taylor_shift(X_Δt, dom)

    # builds the associated taylor model for each coordinate j = 1...n
    #  X̂ is a TaylorModelN whose coefficients are intervals
    X̂ = [TaylorModelN(X_Δt[j], zeroI, zeroBox(n), symBox(n)) for j in 1:n]

    # compute floating point rigorous polynomial approximation
    # fX̂ is a TaylorModelN whose coefficients are floats
    fX̂ = fp_rpa.(X̂)

    # evaluate the spatial variables in the given domain and overapproximate
    Y = _overapproximate_TM(fX̂, T; kwargs...)

    return ReachSet(Y, Δt)
end

# overapproximate a Taylor model with a hyperrectangle
# TODO outsource to LazySets
function _overapproximate_TM(fX̂, ::Type{<:Hyperrectangle}; kwargs...)
    B = IntervalBox([evaluate(fX̂[i], domain(p)) for (i, p) in enumerate(fX̂)]...)
    return convert(Hyperrectangle, B)
end

# overapproximate a Taylor model with a zonotope using LazySets
function _overapproximate_TM(fX̂, ::Type{<:Zonotope}; kwargs...)
    return overapproximate(fX̂, Zonotope; kwargs...)
end

function overapproximate(R::TaylorModelReachSet, T::Type{<:LazySet};
                         Δt::TimeInterval=tspan(R), dom=symBox(dim(R)), kwargs...)
    return _overapproximate(R, T; Δt=Δt, dom=dom, kwargs...)
end

LazySets.box_approximation(R::TaylorModelReachSet) = overapproximate(R, Hyperrectangle)

# overapproximate TaylorModelReachSet with several sets by splitting
# in space (with either nsdiv for uniform partitions or by specifying a
# partition (vector of integers); the version not used should be set to
# `nothing`) and/or in time (ntdiv)
# e.g. if the partition is uniform for each dimension, the number of returned
# sets is nsdiv^D * ntdiv
function overapproximate(R::TaylorModelReachSet{N}, T::Type{<:Union{Hyperrectangle,Zonotope}};
                         partition=nothing, nsdiv=nothing, ntdiv=1,
                         Δt::TimeInterval=tspan(R), dom=symBox(dim(R)), kwargs...) where {N}
    if !isnothing(partition) && !isnothing(nsdiv)
        throw(ArgumentError("either `partition` or `nsdiv` should be specified, not both"))
    end
    if isnothing(partition) && isnothing(nsdiv)
        nsdiv = 1
    end

    # no splitting
    if isnothing(partition) && nsdiv == 1 && ntdiv == 1
        return _overapproximate(R, T; Δt=Δt, dom=dom, kwargs...)
    end

    D = dim(R)
    S = symBox(D)
    @assert Δt == tspan(R) "time subdomain not implemented"
    @assert dom == S "spatial subdomain not implemented"

    X = set(R)

    # time domain splittig
    tspdiv = IA.mince(tspan(R), ntdiv)
    tdom = domain(R)
    Δtdiv = IA.mince(tdom, ntdiv)

    # evaluate the Taylor model in (normalized) time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = [evaluate(X, Δti) for Δti in Δtdiv]

    # partition the symmetric box
    if isnothing(partition)
        if nsdiv == 1
            partition = [S]
        else
            # TODO (also below) may use IA.mince directly
            partition = fill(nsdiv, D)
            Sdiv = LazySets.split(convert(T, S), partition)
            partition = convert.(IntervalBox, Sdiv)
        end

    else
        Sdiv = LazySets.split(convert(T, S), partition)
        partition = convert.(IntervalBox, Sdiv)
    end
    nparts = length(partition)

    # preallocate output vectors
    if T <: Hyperrectangle
        ST = Hyperrectangle{N,SVector{D,N},SVector{D,N}}
        R̂out = Vector{ReachSet{N,ST}}(undef, nparts * ntdiv)
    else
        R̂out = Vector{ReachSet{N}}(undef, nparts * ntdiv)
    end

    # evaluate the spatial variables in the symmetric box
    @inbounds for k in 1:ntdiv
        for j in 1:nparts
            Bj = partition[j]
            Xk = X_Δt[k]
            X̂ib = IntervalBox([evaluate(Xk[i], Bj) for i in 1:D])

            idx = j + (k - 1) * nparts
            X̂out = convert(T, convert(Hyperrectangle, X̂ib))

            R̂out[idx] = ReachSet(X̂out, tspdiv[k])
        end
    end
    if T <: Zonotope  # concrete output vector
        R̂out = [e for e in R̂out]
    end

    return R̂out
end

# overapproximate taylor model reachset with several zonotopes given the partition,
# that specifies the way in which the symmetric box [-1, 1]^n is split:
#
# - if `partition` is an integer, makes as uniform partition for each coordinate
# - if `partition` is a vector of integers, split the domain according to each element in partition
function overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope},
                         partition::Union{Int,Vector{Int}}) where {N}
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
    fX̂ = Vector{Vector{TaylorModelN{length(X_Δt),N,N}}}(undef, length(part))
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

# evaluate at a given time point and overapproximate the resulting set
function overapproximate(R::TaylorModelReachSet, T::Type{<:LazySet}, t::AbstractFloat; kwargs...)
    return overapproximate(R, T; Δt=TimeInterval(t, t), kwargs...)
end

# convert a hyperrectangular set to a taylor model reachset
function convert(::Type{<:TaylorModelReachSet}, H::AbstractHyperrectangle{N};
                 orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI) where {N}
    n = dim(H)
    x = set_variables("x"; numvars=n, order=orderQ)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N},N}}(undef, n)

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
    return convert(T, H; orderQ=orderQ, orderT=orderT, Δt=Δt)
end

# overapproximate a zonotopic set with a taylor model reachset
# TODO pass algorithm option to choose between using parallelotope oa or order reduction
function overapproximate(Z::AbstractZonotope{N}, ::Type{<:TaylorModelReachSet};
                         orderQ::Integer=2, orderT::Integer=8, Δt::TimeInterval=zeroI,
                         indices=1:dim(Z), box_reduction=false) where {N}
    n = dim(Z)
    x = set_variables("x"; numvars=n, order=orderQ)

    if order(Z) > 1
        if box_reduction
            # diagonal generators matrix
            Z = reduce_order(Z, 1)
        else
            # indices selects the indices that we want to keep
            Z = _overapproximate_hparallelotope(Z, indices)
        end
    end
    c = center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N},N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # for each variable i = 1, .., n, compute the linear polynomial that covers
    # the line segment corresponding to the i-th edge of Z
    p = size(G, 2)  # the order may be < 1
    @inbounds for i in 1:n
        pi = c[i] + sum(view(G, i, :) .* x[1:p])
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
                                     orderQ::Integer=2, orderT::Integer=8,
                                     Δt::TimeInterval=zeroI) where {N}
    n = dim(Z)
    x = set_variables("x"; numvars=n, order=orderQ)

    # check structure
    order(Z) == 2 ||
        throw(ArgumentError("this function requires that the order of the zonotope is 2, got $(order(Z))"))

    c = LazySets.center(Z)
    G = genmat(Z)

    M = view(G, :, 1:n)
    D = view(G, :, (n + 1):(2n))
    isdiag(D) ||
        throw(ArgumentError("the columns $(n+1) to $(2n) of the generators matrix do not form a diagonal matrix"))

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N},N}}(undef, n)

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

function _overapproximate_structured(Zcp::CartesianProduct{N,<:Zonotope,<:Interval},
                                     ::Type{<:TaylorModelReachSet};
                                     orderQ::Integer=2, orderT::Integer=8,
                                     Δt::TimeInterval=zeroI) where {N}
    n = dim(Zcp)
    x = set_variables("x"; numvars=n, order=orderQ)

    # check structure
    Z = Zcp.X
    Y = Zcp.Y
    if order(Z) > 1
        Z = remove_redundant_generators(Z)
    end
    order(Z) == 1 ||
        throw(ArgumentError("this function requires that the order of the (lower) zonotope is 1, got $(order(Z))"))

    c = LazySets.center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N},N}}(undef, n)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    xstate = x[1:(n - 1)]
    @inbounds for i in 1:(n - 1)
        pi = c[i] + sum(view(G, i, :) .* xstate) + zero(TaylorN(n; order=orderQ))
        rem = interval(0)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end
    # diagonal part
    @inbounds begin
        pi = mid(Y.dat) + zero(TaylorN(1; order=orderQ))
        d = diam(Y.dat) / 2
        rem = interval(-d, d)
        vTM[n] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end
    return TaylorModelReachSet(vTM, Δt)
end

# we assume that if n = size(G, 1), then:
# (size(G) == (n, 2n + 1)) && isdiag(view(G, :, (n+2):(2n+1)))
function _overapproximate_structured_full(Zcp::CartesianProduct{N,<:Zonotope,<:Interval},
                                          ::Type{<:TaylorModelReachSet};
                                          orderQ::Integer=2, orderT::Integer=8,
                                          Δt::TimeInterval=zeroI) where {N}
    n = dim(Zcp) - 1
    x = set_variables("x"; numvars=n + 1, order=orderQ)

    # check structure
    # not checking structure
    Z = Zcp.X
    c = LazySets.center(Z)
    G = genmat(Z)

    # preallocations
    vTM = Vector{TaylorModel1{TaylorN{N},N}}(undef, n + 1)

    # normalized time domain
    Δtn = TimeInterval(zero(N), diam(Δt))

    # fill rows corresponding to the "zonotope" variables: 1 to nth-variables
    @inbounds for i in 1:n
        pi = c[i] + sum(view(G, i, 1:(n + 1)) .* x) + zero(TaylorN(n + 1; order=orderQ))
        d = abs(G[i, n + 1 + i])
        rem = interval(-d, d)
        vTM[i] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)
    end

    # fill the final row, which correspponds to the "interval" variable: (n+1)-th
    I = Zcp.Y.dat
    pi = mid(I) + zero(TaylorN(n + 1; order=orderQ))
    d = diam(I) / 2
    rem = interval(-d, d)
    @inbounds vTM[n + 1] = TaylorModel1(Taylor1(pi, orderT), rem, zeroI, Δtn)

    return TaylorModelReachSet(vTM, Δt)
end
