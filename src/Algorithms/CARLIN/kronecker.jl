# ======================================
# Kronecker powers
# ======================================

using LinearAlgebra: checksquare

"""
    kron_pow(x::IA.Interval, pow::Int)

Given an interval x and an integer pow, compute `x^pow`.

### Input

- `x`   -- interval
- `pow` -- integer

### Output

An interval enclosure of `x^pow`.

### Examples

```julia
julia> [kron_pow(2 .. 3, i) for i in 1:3]
3-element Array{IntervalArithmetic.Interval{Float64},1}:
  [2, 3]
  [4, 9]
 [8, 27]
```
"""
function kron_pow(x::IA.Interval, pow::Int)
    return x^pow
end

# compute x^{⊗ pow} where x is a LazySets.Interval
function kron_pow(x::Interval, pow::Int)
    return Interval(kron_pow(x.dat, pow))
end

# compute H^{⊗ pow} where H is a hyperrectangular set
# we see H as a product of intervals and apply
# kron_pow(x, pow) for each interval element x of H
function kron_pow(H::AbstractHyperrectangle, pow::Int)
    Hbox = convert(IntervalBox, H)
    Hbox_pow = kron_pow.(Hbox, pow)
    #Hbox_pow = convert(Hyperrectangle, Hbox_pow)
    #Hout = Hyperrectangle(Vector(Hbox_pow.center), Vector(Hbox_pow.radius))
    return convert(Hyperrectangle, Hbox_pow)
end

# compute [x, x^2, .., x^pow] where x is an interval
function kron_pow_stack(x::Union{<:Interval, <:IA.Interval}, pow::Int)
    return [kron_pow(x, i) for i in 1:pow]
end

# compute the cartesian product [x, x^2, .., x^pow] where x is a hyperrectangular set
function kron_pow_stack(H::AbstractHyperrectangle, pow::Int)
    out = [kron_pow(H, i) for i in 1:pow]
    return CartesianProductArray(out)
end

"""
    kron_id(n, k)

Compute the k-fold Kronecker product of the identity: `I ⊗ I ⊗ ... ⊗ I`, k times

### Input

- `n` -- integer representing the order (dimension of the identity)
- `k` -- integer representing the power

### Output

A `Diagonal` matrix with element `1` in the diagonal and order ``n^k``.

### Examples

```julia
julia> kron_id(2, 2)
4×4 Diagonal{Float64,Array{Float64,1}}:
 1.0   ⋅    ⋅    ⋅
  ⋅   1.0   ⋅    ⋅
  ⋅    ⋅   1.0   ⋅
  ⋅    ⋅    ⋅   1.0
```
"""
function kron_id(n::Int, k::Int)
    return Diagonal(ones(n^k))
end

"""
    kron_sandwich(F, n, k1, k2)

Compute `A ⊗ F ⊗ C` where `A = I^{⊗ k1}` and `C = I^{⊗ k2}` where `I` is the
identity matrix of order `n`, the `k1`-fold and `k2`-fold Kronecker product of the
identity and `F` in between.

### Input

- `F`  -- matrix
- `n`  -- integer, dimension of the identity
- `k1` -- nonnegative integer
- `k2` -- nonnegative integer

### Output

The kronecker product `I^{⊗ k1} ⊗ F ⊗ I^{⊗ k2}`, represented as a sparse matrix
if `F` is sparse.

### Examples

```julia
julia> F = sparse([0 1; -1 0.])
2×2 SparseMatrixCSC{Float64,Int64} with 2 stored entries:
  [2, 1]  =  -1.0
  [1, 2]  =  1.0

julia> Q = kron_sandwich(F, 2, 2, 2);

julia> size(Q)
(32, 32)

julia> I2 = I(2)
IdentityMultiple{Float64} of value 1.0 and order 2

julia> Q == reduce(kron, [I2, I2, F, I2, I2])
true
```
"""
function kron_sandwich(F::AbstractMatrix, n::Int, k1::Int, k2::Int)
    # compute A ⊗ F
    if k1 == 0
        AF = F
    else
        A = kron_id(n, k1)
        AF = kron(A, F)
    end

    # compute A ⊗ F ⊗ C
    if k2 == 0
        AFC = AF
    else
        C = kron_id(n, k2)
        AFC = kron(AF, C)
    end

    return AFC
end

"""
    kron_sum(F::AbstractMatrix, k::Int)

Compute the Kronecker sum of order k defined as:

```math
    F ⊕ F ⊕ ... ⊕ F := F ⊗ I ⊗ ... ⊗ I + I ⊗ F ⊗ I ⊗ ... ⊗ I + ... + I ⊗ ... ⊗ I ⊗ F
```
where each term has `k` products and there are a total of `k` summands
and `I` is the identity matrix of order `n`.

### Examples:

It holds that:

- `k = 1: F`
- `k = 2: F ⊗ I + I ⊗ F`
- `k = 3: F ⊗ I ⊗ I + I ⊗ F ⊗ I + I ⊗ I ⊗ F`

```julia
julia> F = sparse([0 1; -1 0.])
2×2 SparseMatrixCSC{Float64,Int64} with 2 stored entries:
  [2, 1]  =  -1.0
  [1, 2]  =  1.0

julia> kron_sum(F, 1) == F
true

julia> I2 = I(2)
IdentityMultiple{Float64} of value 1.0 and order 2

julia> kron_sum(F, 2) == kron(F, I2) + kron(I2, F)
true

julia> kron_sum(F, 3) == sum(reduce(kron, X) for X in [[F, I2, I2], [I2, F, I2], [I2, I2, F]])
true
"""
function kron_sum(F::AbstractMatrix, k::Int)
    if k < 1
        throw(ArgumentError("expected k ≥ 1, got $k"))
    elseif k == 1
        return F
    end

    n = checksquare(F)

    k1 = 0 # terms on the left of F
    k2 = k-1 # terms on the right of F
    A = kron_sandwich(F, n, k1, k2)  # I^(⊗k1) ⊗ F ⊗ I^(⊗k2)
    for i in 2:k
        k1 += 1
        k2 -= 1
        B = kron_sandwich(F, n, k1, k2)
        A += B
    end
    return A
end

# compute matrix A =
# ------------------------------------------------------------
# Functionality that requires MultivariatePolynomials.jl
# ------------------------------------------------------------

function load_kron_multivariate()
return quote

import Base: findfirst, findall

"""
    kron_pow(x::Vector{<:AbstractVariable}, pow::Int)

Compute the higher order concrete Kronecker power: `x ⊗ x ⊗ ... ⊗ x`, `pow` times
for a vector of symbolic monomials.

### Input

- `x`   -- polynomial variable
- `pow` -- integer

### Output

Vector of multivariate monomial corresponding to `x^{⊗ pow}`.

### Examples

```julia
julia> using DynamicPolynomials

julia> @polyvar x[1:2]
(PolyVar{true}[x₁, x₂],)

julia> x
2-element Array{PolyVar{true},1}:
 x₁
 x₂

julia> kron_pow(x, 2)
4-element Array{Monomial{true},1}:
 x₁²
 x₁x₂
 x₁x₂
 x₂²
 ```
"""
function kron_pow(x::Vector{<:AbstractVariable}, pow::Int)
    @assert pow > 0 "expected positive power, got $pow"
    if pow == 1
        return x
    else
        return kron(x, kron_pow(x, pow-1))
    end
end

"""
    findfirst(y::Vector{<:AbstractMonomialLike}, x::AbstractMonomialLike)

Return the first position of the multivariate monomial ``x`` in the vector of monomials ``y``.

### Input

- `y` -- vector of multivariate monomials
- `x` -- multivariate monomials

### Output

An integer where each integer represents the index in the array `y`
corresponding to the first match, i.e. `y[i] == x` where `i` is the output integer,
or `nothing` if there is no match.

### Notes

Let ``x = (x_1, x_2, …, x_n)`` be given, and let ``y_i = x^{[i]}``.
Given the multi-index ``I = (i_1, i_2, …, i_n)``, this function returns the
first position of ``x^I`` in the array ``y`` (resp. all positions using `findall`).

### Examples

```julia
julia> using DynamicPolynomials

julia> @polyvar x[1:2]
(PolyVar{true}[x₁, x₂],)

julia> y = kron_pow(x, 2)
4-element Array{Monomial{true},1}:
 x₁²
 x₁x₂
 x₁x₂
 x₂²

julia> findfirst(y, x[1]*x[2])
2
```
"""
function findfirst(y::Vector{<:AbstractMonomialLike}, x::AbstractMonomialLike)
    pow = exponents(x)
    @assert sum(pow) == sum(exponents(first(y))) "power indices don't match the power in the lifted vector"
    xv = variables(x)
    _findfirst(y, xv, pow)
end

"""
    findall(y::Vector{<:AbstractMonomialLike}, x::AbstractMonomialLike)

Return all positions of the multivariate monomial ``x`` in the vector of monomials ``y``.

### Input

- `y` -- vector of multivariate monomials
- `x` -- multivariate monomials

### Output

A vector of integers where each integer represents the index in the array `y`
corresponding to a match, i.e. `y[i] == x` for all elements `i` in the output array.

### Notes

Let ``x = (x_1, x_2, …, x_n)`` be given, and let ``y_i = x^{[i]}``.
Given the multi-index ``I = (i_1, i_2, …, i_n)``, this function returns all
positions of ``x^I`` in the array ``y`` (resp. the first position using `findfirst`).

### Examples

```julia
julia> using DynamicPolynomials

julia> @polyvar x[1:2]
(PolyVar{true}[x₁, x₂],)

julia> y = kron_pow(x, 2)
4-element Array{Monomial{true},1}:
 x₁²
 x₁x₂
 x₁x₂
 x₂²

julia> findall(y, x[1]*x[2])
2-element Array{Int64,1}:
 2
 3
```
"""
function findall(y::Vector{<:AbstractMonomialLike}, x::AbstractMonomialLike)
    pow = exponents(x)
    @assert sum(pow) == sum(exponents(first(y))) "power indices don't match the power in the lifted vector"
    xv = variables(x)
    _findall(y, xv, pow)
end

function _findfirst(y, x, pow)
    yp = powers.(y)
    t = ntuple(i -> (x[i], pow[i]), length(pow))
    idx = findfirst(pi -> pi == t, yp)
end

function _findall(y, x, pow)
    yp = powers.(y)
    t = ntuple(i -> (x[i], pow[i]), length(pow))
    idx = findall(pi -> pi == t, yp)
end

end end  # quote / load_kron_multivariate()
