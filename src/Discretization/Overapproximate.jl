module Overapproximate

using LazySets
using IntervalMatrices
using LazySets: matrix
using StaticArrays: SVector, SMatrix, MMatrix, StaticArray
import IntervalArithmetic as IA

export _convert_or_overapproximate, _overapproximate, _split, relative_error,
       _overapproximate_hparallelotope

# ===========================
# _convert_or_overapproximate
# ===========================

# fallback implementation for conversion (if applicable) or overapproximation
function _convert_or_overapproximate(T::Type{<:AbstractPolytope}, X::LazySet)
    if applicable(convert, T, X)
        return convert(T, X)
    end
    try
        return overapproximate(X, T)
    catch e
        if !(e isa MethodError)
            rethrow()
        end
    end
    return convert(T, overapproximate(X, Hyperrectangle))
end

# ================
# _overapproximate
# ================

function _overapproximate(X::AbstractPolytope, T::Type{<:HPolytope})
    P = HPolytope(constraints_list(X))
    return convert(T, P)
end

# overapproximate a hyperrectangular set with a polytope with `Vector` directions
# (by default it uses `SingleEntryVector` directions)
function _overapproximate(H::AbstractHyperrectangle, ::Type{<:HPolytope})
    return HPolytope([HalfSpace(Vector(c.a), c.b) for c in constraints_list(H)])
end

function _overapproximate(lm::LinearMap{N,<:AbstractZonotope{N},NM,<:AbstractIntervalMatrix{NM}},
                          ::Type{<:Zonotope}) where {N<:Real,NM}
    Mc, Ms = _split(matrix(lm))
    Z = LazySets.set(lm)
    c = center(Z)
    G = genmat(Z)
    return _overapproximate_interval_linear_map(Mc, Ms, c, G)
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
    q = size(Ms, 1)
    α = Ms * dvec # vector of length q
    αnz = findall(!iszero, α)
    DV = zeros(N, q, length(αnz))
    @inbounds for (j, idx) in enumerate(αnz)
        DV[idx, j] = α[idx]
    end
    G_oa = hcat(Ggens, DV)
    return Zonotope(c_oa, G_oa)
end

function _overapproximate_interval_linear_map(Mc::SMatrix{n,n,N,LM},
                                              Ms::SMatrix{n,n,N,LM},
                                              c::SVector{n,N},
                                              G::SMatrix{n,m,N,LG}) where {n,N,LM,m,LG}
    c_oa = Mc * c
    Ggens = Mc * G

    dvec = zeros(N, n)
    @inbounds for i in 1:n
        dvec[i] = abs(c[i])
        for j in 1:m
            dvec[i] += abs(G[i, j])
        end
    end
    q = size(Ms, 1)
    α = Ms * dvec # vector of length q
    αnz = findall(!iszero, α)
    DV = zeros(MMatrix{q,q,N}) # NOTE: sole difference with regular arrays, may refactor
    @inbounds for (j, idx) in enumerate(αnz)
        DV[idx, j] = α[idx]
    end
    G_oa = hcat(Ggens, DV)
    return Zonotope(c_oa, G_oa)
end

# type-stable version
function _overapproximate(S::LazySet{N}, ::Type{<:Hyperrectangle}) where {N}
    c, r = box_approximation_helper(S)
    #if r[1] < 0
    #    return EmptySet{N}(dim(S))
    #end
    return Hyperrectangle(c, r)
end

# compared to LazySets.Approximations._overapproximate_hparallelotope,
# this function does inv(Matrix(Γ))
function _overapproximate_hparallelotope(Z::AbstractZonotope, indices=1:dim(Z))
    length(indices) == dim(Z) ||
        throw(ArgumentError("the number of generator indices is $(length(indices)), " *
                            "but it was expected to be $(dim(Z))"))

    p, n = ngens(Z), dim(Z)
    if p == n
        return Z
    elseif p < n
        error("the zonotope order is $(order(Z)) but it should be at least 1")
    end

    G = genmat(Z)
    Γ = G[:, indices]
    □Γ⁻¹Z = box_approximation(linear_map(inv(Matrix(Γ)), Z))
    return linear_map(Γ, □Γ⁻¹Z)
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

function _split(A::IntervalMatrix{T,IT,MT}) where {T,IT,MT<:AbstractMatrix{IT}}
    m, n = size(A)
    C = Matrix{T}(undef, m, n)
    S = Matrix{T}(undef, m, n)
    _split_fallback!(A, C, S)
    return C, S
end

function _split(A::AbstractMatrix)
    return A, zeros(size(A))
end

function _split(A::IntervalMatrix{T,IT,MT}) where {T,IT,ST,MT<:StaticArray{ST,IT}}
    m, n = size(A)
    # TODO: use MMatrix and convert to SMatrix afterwards?
    C = Matrix{T}(undef, m, n)
    S = Matrix{T}(undef, m, n)
    _split_fallback!(A, C, S)
    return SMatrix{m,n,T}(C), SMatrix{m,n,T}(S)
end

# TEMP
function LazySets.Approximations.box_approximation(x::IA.Interval)
    return convert(Hyperrectangle, Interval(x))
end

# TEMP
function LazySets.Approximations.box_approximation(x::IA.IntervalBox)
    return convert(Hyperrectangle, x)
end

@inline function box_approximation_helper(S::LazySet{N}) where {N}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        l, h = extrema(S, i)
        r[i] = (h - l) / 2
        if r[i] < 0
            # contradicting bounds => set is empty
            # terminate with first radius entry being negative
            r[1] = r[i]
            break
        end
        c[i] = (h + l) / 2
    end
    return c, r
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
(see [AlthoffGK18; Eq. (20)](@citet)).
"""
function relative_error(x, x_ref)
    x_low, x_high = inf(x), sup(x)
    x_ref_low, x_ref_high = inf(x_ref), sup(x_ref)
    denom = x_ref_high - x_ref_low
    rel_low = -(x_low - x_ref_low) / denom
    rel_high = (x_high - x_ref_high) / denom
    return 100 * IA.interval(rel_low, rel_high)
end

end  # module
