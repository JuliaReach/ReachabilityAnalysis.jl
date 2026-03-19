"""
Interface to matrix exponential backends of different kinds.

Includes common integral computations arising in the discretization of linear
differential equations using matrix methods. For applications see, e.g.,
[ForetsS22](@citet) and references therein.
"""
module Exponentiation

using StaticArrays: StaticArray
using MathematicalSystems: IdentityMultiple
using ReachabilityBase.Require: require
using LinearAlgebra: checksquare, I
using SparseArrays: SparseMatrixCSC, nonzeros, sparse, spzeros
using IntervalMatrices: AbstractIntervalMatrix, IntervalMatrix,
                        exp_overapproximation
using LazySets: Hyperrectangle, LazySet, SparseMatrixExp, dim,
                get_columns, symmetric_interval_hull, _expmv

export BaseExp, BaseExpAlg, IntervalExpAlg, LazyExpAlg, PadeExpAlg, elementwise_abs, Φ₂, Φ₁, Φ₁_u

# -------------------------
# Exponentiation interface

"""
    AbstractExpAlg

Abstract supertype for all exponentiation algorithms.
"""
abstract type AbstractExpAlg end

# TODO use _alias_exp to distinguish from setops
# no-op
_alias(alg::AbstractExpAlg) = alg

# alias defined on value types
_alias(alg::Symbol) = _alias(Val(alg))

"""
    BaseExpAlg <: AbstractExpAlg

Matrix exponential using the scaling-and-squaring algorithm implemented in Julia
Base.

### Notes

The alias for this algorithm is `:base`.
"""
struct BaseExpAlg <: AbstractExpAlg end

const BaseExp = BaseExpAlg()
_alias(::Val{:base}) = BaseExp

"""
    LazyExpAlg <: AbstractExpAlg

Matrix exponential computed in a lazy fashion.

### Fields

- `m`   -- (optional, default: `30`) size of the Krylov subspace
- `tol` -- (optional, default: `1e-10`) tolerance

### Notes

The aliases for this algorithm are `:lazy` and `:krylov`.

The operations are defined in the package `LazySets.jl` (`SparseMatrixExp`).
"""
struct LazyExpAlg <: AbstractExpAlg
    m::Int
    tol::Float64
end

# default constructor
LazyExpAlg() = LazyExpAlg(30, 1e-10)

const LazyExp = LazyExpAlg()
_alias(::Union{Val{:lazy},Val{:krylov}}) = LazyExp

"""
    IntervalExpAlg <: AbstractExpAlg

Matrix exponential using an interval enclosure of the Taylor series remainder.

### Fields

- `order` -- order of the Taylor series expansion of the matrix exponential

### Notes

The aliases for this algorithm are `:interval` and `:taylor`.

This algorithm allows to overapproximate ``exp(Aδ)`` with an interval matrix.
It also accepts an interval matrix ``A`` as input.
"""
struct IntervalExpAlg <: AbstractExpAlg
    order::Int
end

# default constructor
IntervalExpAlg() = IntervalExpAlg(10)

const IntervalExp = IntervalExpAlg()
_alias(::Union{Val{:interval},Val{:taylor}}) = IntervalExp

"""
    PadeExpAlg <: AbstractExpAlg

Matrix exponential for sparse matrices using Pade approximants.

### Notes

The alias for this algorithm is `:pade`.

This algorithm requires to load the package `Expokit.jl`.
"""
struct PadeExpAlg <: AbstractExpAlg end

const PadeExp = PadeExpAlg()
_alias(::Val{:pade}) = PadeExp

"""
    _exp(A::AbstractMatrix, δ, [alg]::AbstractExpAlg=BaseExp)

Compute the matrix exponential ``e^{Aδ}``.

### Input

- `A`    -- matrix
- `δ`    -- step size
- `alg`  -- (optional, default: `BaseExp`) the algorithm used to take the matrix
            exponential of `Aδ`, possible options are `BaseExp`, `LazyExp`,
            `PadeExp` and `IntervalExp` (see details in the *Algorithm* section
            below)

### Output

A matrix or a lazy wrapper of the matrix exponential, depending on `alg`.

### Algorithm

- `BaseExp` -- (alias: `:base`) use Higham's scaling-and-squaring method implemented
               in Julia's standard library; see `?exp` for details; if `A` is a static array,
               uses the implementation in `StaticArrays.jl`

- `LazyExp` -- (alias: `:lazy`) return a lazy wrapper around the matrix exponential
               using the implementation `LazySets.SparseMatrixExp`

- `PadeExp` -- (alias: `pade`) apply the Padé approximant method to compute the matrix
               exponential of a sparse matrix (requires `Expokit.jl`)

- `IntervalExp` -- (alias: `interval`, `taylor`) apply the Taylor series expansion of the matrix
                   exponential with an interval remainder; works if `A` is an interval matrix

### Notes

If the algorithm `LazyExp` is used, actions of the matrix exponential are
evaluated with an external library such as `ExponentialUtilities.jl` or
`Expokit.jl`.
"""
function _exp(A::AbstractMatrix, δ, alg::AbstractExpAlg=BaseExp)
    checksquare(A)
    return _exp(A * δ, alg)
end

# general case: convert to Matrix
_exp(A::AbstractMatrix, ::BaseExpAlg) = exp(Matrix(A))

# Base's `exp`
_exp(A::Matrix, ::BaseExpAlg) = exp(A)

# static arrays have their own exp method
_exp(A::StaticArray, ::BaseExpAlg) = exp(A)

# exponential of an identity multiple, defined in MathematicalSystems.jl
_exp(A::IdentityMultiple, ::BaseExpAlg) = exp(A)

# lazy wrapper of the matrix exponential
_exp(A::AbstractMatrix, ::LazyExpAlg) = SparseMatrixExp(A)

# pade approximants (requires Expokit.jl)
_exp(::AbstractMatrix, ::PadeExpAlg) = throw(ArgumentError("algorithm requires a sparse matrix"))

function _exp(A::SparseMatrixCSC, ::PadeExpAlg)
    require(@__MODULE__, :Expokit)
    return _exp_pade(A)
end

function load_expokit_pade()
    return quote
        _exp_pade(A::SparseMatrixCSC) = Expokit.padm(A)
    end
end  # quote / load_expokit_pade

function _exp(A::AbstractIntervalMatrix, δ, alg::IntervalExpAlg)
    return exp_overapproximation(A, δ, alg.order)
end

# convert to IntervalMatrix
_exp(A::AbstractMatrix, δ, alg::IntervalExpAlg) = _exp(IntervalMatrix(A), δ, alg)

# add δ explicitly
_exp(A::AbstractMatrix, alg::IntervalExpAlg) = _exp(A, one(eltype(A)), alg)

# ------------------------
# Discretization functions
# ------------------------

"""
    Φ₁(A::AbstractMatrix, δ, [alg]::AbstractExpAlg=BaseExp, [isinv]::Bool=false, [Φ]=nothing)

Evaluate the series

```math
Φ₁(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+1}}{(i+1)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ > 0``.

### Input

- `A`      -- coefficients matrix
- `δ`      -- step size
- `alg`    -- (optional, default: `BaseExp`) the method used to take the matrix
               exponential of the coefficient matrix; see the documentation of
               [`_exp`](@ref) for available options
- `isinv`  -- (optional, default: `false`) if `true`, assume that the coefficients
              matrix is invertible and compute ``A^{-1}``
- `Φ`      -- (optional, default: `nothing`) optionally pass the matrix exponential ``e^{Aδ}``

### Output

A matrix.

### Algorithm

If ``A`` is invertible, ``Φ₁`` can be computed as

```math
Φ₁(A, δ) = A^{-1}(e^{Aδ} - I_n),
```
where ``I_n`` is the identity matrix of order ``n``.

In the general case, implemented in this function, it can be
computed as submatrices of the block matrix

```math
P_{2n} = \\exp \\begin{pmatrix}
Aδ && δI_n \\\\
0 && 0
\\end{pmatrix}.
```
It can be shown that
```math
\\exp(P_{2n}) = \\begin{pmatrix}
Φ(A, δ) && Φ₁(A, δ) \\\\
0 &&  δI_n
\\end{pmatrix}.
```
where ``Φ(A, δ) = e^{Aδ}``. In particular, `Φ₁(A, δ) = P[1:n, (n+1):2*n]`.
This method can be found in [FrehseGDCRLRGDM11](@citet).
"""
function Φ₁(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg=BaseExp, isinv::Bool=false, Φ=nothing)
    return Φ₁(A, δ, alg, Val(isinv), Φ)
end

# dispatch for discretization methods
Φ₁(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg, ::Val{:false}, Φ) = _Φ₁_blk(A, δ, alg)
Φ₁(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg, ::Val{:true}, Φ) = _Φ₁_inv(A, δ, alg, Φ)

# evaluate the series Φ₁(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+1}}{(i+1)!}A^i
# without assuming invertibility of A, by taking the exponential of a 2n x 2n matrix
function _Φ₁_blk(A, δ, alg)
    n = checksquare(A)
    P_2n = _P_2n(A, δ, n)
    Q = _exp(P_2n, alg)
    return _P₁_blk(Q, n, alg)
end

@inline function _P_2n(A::AbstractMatrix{N}, δ, n) where {N}
    return [Matrix(A * δ) Matrix(δ * I, n, n);
            zeros(n, 2 * n)]::Matrix{N}
end

@inline function _P_2n(A::SparseMatrixCSC{N,M}, δ, n) where {N,M}
    return [sparse(A * δ) sparse(δ * I, n, n);
            spzeros(n, 2 * n)]::SparseMatrixCSC{N,M}
end

@inline _P₁_blk(P, n::Int, ::AbstractExpAlg) = P[1:n, (n + 1):(2 * n)]
@inline _P₁_blk(P, n::Int, ::LazyExpAlg) = sparse(get_columns(P, (n + 1):(2 * n))[1:n, :])

# compute the matrix Φ₁(A, δ) = A^{-1}(exp(Aδ) - I), assuming that A is invertible
# and explicitly computing inv(A); this function optionally receives the matrix Φ = exp(Aδ)
function _Φ₁_inv(A::AbstractMatrix, δ, alg, Φ=nothing)
    if isnothing(Φ)
        Φ = _exp(A, δ, alg)
    end
    n = size(A, 1)
    N = eltype(A)
    In = Matrix(one(N) * I, n, n)
    Ainv = inv(A)
    return Ainv * (Φ - In)
end

# special case when A is a multiple of the identity
function _Φ₁_inv(A::IdentityMultiple, δ, alg, Φ=nothing)
    λ = A.M.λ
    @assert !iszero(λ) "the given identity multiple is not invertible"
    α = (1 / λ) * (exp(δ * λ) - 1)
    return IdentityMultiple(α, size(A, 1))
end

# method to compute Φ₁ * u, where u is a vector
function Φ₁_u(A, δ, alg, ::Val{:true}, u::AbstractVector, Φ=nothing)
    return Φ₁_u(A, δ, alg, Φ, u)
end

function Φ₁_u(A, δ, alg, ::Val{:false}, u::AbstractVector, Φ=nothing)
    M = _Φ₁_blk(A, δ, alg)
    return M * u
end

# compute the vector Φ₁(A, δ)u = A^{-1}(exp(Aδ) - I) u
# assuming that A is invertible and solving the associated linear system
# (the inverse of A is not explicitly computed)
function Φ₁_u(A, δ, alg, u::AbstractVector, Φ=nothing)
    if isnothing(Φ)
        Φ = _exp(A, δ, alg)
    end
    x = Φ * u - u
    return A \ x
end

# compute Φ₁(A, δ)u = A^{-1}(exp(Aδ) - I) u without explicitly computing exp(Aδ)
# and assuming that A is invertible
function Φ₁_u(A, δ, alg::LazyExpAlg, u::AbstractVector, ::Nothing)
    w = _expmv(δ, A, u; m=alg.m, tol=alg.tol)
    x = w - u
    return A \ x
end

"""
    Φ₂(A::AbstractMatrix, δ, [alg]::AbstractExpAlg=BaseExp, [isinv]::Bool=false, [Φ]=nothing)

Evaluate the series

```math
Φ₂(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+2}}{(i+2)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ > 0``.

### Input

- `A`      -- coefficients matrix
- `δ`      -- step size
- `alg`    -- (optional, default: `BaseExp`) the method used to take the matrix
               exponential of the coefficient matrix; see the documentation of
               [`_exp`](@ref) for available options
- `isinv`  -- (optional, default: `false`) if `true`, assume that the coefficients
              matrix is invertible and compute ``A^{-1}``
- `Φ`      -- (optional, default: `nothing`) optionally pass the matrix exponential ``e^{Aδ}``

### Output

A matrix.

### Algorithm

If ``A`` is invertible, ``Φ₂`` can be computed as

```math
Φ₂(A, δ) = A^{-2}(e^{δA} - I_n - δA).
```

In the general case, implemented in this function, it can be computed as
submatrices of the block matrix

```math
P_{3n} = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\\\
0 && 0 && δI_n \\\\
0 && 0 && 0
\\end{pmatrix}.
```
It can be shown that
```math
\\exp(P_{3n}) = \\begin{pmatrix}
Φ(A, δ) && Φ₁(A, δ) && Φ₂(A, δ) \\\\
0 && I_n     && δI_n \\\\
0 && 0     && I_n
\\end{pmatrix}.
```
where ``Φ(A, δ) = e^{Aδ}``. In particular, `Φ₂ = P_{3n}[1:n, (2*n+1):3*n]`.
This method can be found in [FrehseGDCRLRGDM11](@citet).
"""
function Φ₂(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg=BaseExp, isinv::Bool=false, Φ=nothing)
    return Φ₂(A, δ, alg, Val(isinv), Φ)
end

# dispatch for discretization methods
Φ₂(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg, ::Val{:false}, Φ) = _Φ₂_blk(A, δ, alg)
Φ₂(A::AbstractMatrix, δ::Real, alg::AbstractExpAlg, ::Val{:true}, Φ) = _Φ₂_inv(A, δ, alg, Φ)

# evaluate the series Φ₂(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+2}}{(i+2)!}A^i
# without assuming invertibility of A, by taking the exponential of a 3n x 3n matrix
function _Φ₂_blk(A, δ, alg)
    n = checksquare(A)
    B = _P_3n(A, δ, n)
    P = _exp(B, alg)
    return _P₂_blk(P, n, alg)
end

@inline function _P_3n(A::AbstractMatrix{N}, δ, n) where {N}
    return [Matrix(A * δ) Matrix(δ * I, n, n) zeros(n, n);
            zeros(n, 2 * n) Matrix(δ * I, n, n);
            zeros(n, 3 * n)]::Matrix{N}
end

@inline function _P_3n(A::SparseMatrixCSC{N,M}, δ, n) where {N,M}
    return [sparse(A * δ) sparse(δ * I, n, n) spzeros(n, n);
            spzeros(n, 2 * n) sparse(δ * I, n, n);
            spzeros(n, 3 * n)]::SparseMatrixCSC{N,M}
end

@inline _P₂_blk(P, n, ::AbstractExpAlg) = P[1:n, (2 * n + 1):(3 * n)]
@inline _P₂_blk(P, n, ::LazyExpAlg) = sparse(get_columns(P, (2 * n + 1):(3 * n))[1:n, :])

# compute the matrix Φ₂ = A^{-2} (exp(A*δ) - I - A*δ) assuming that A is invertible
# and explicitly computing inv(A); this function optionally receives Φ = exp(Aδ)
# Φ2 = ReachabilityAnalysis._Φ₂_inv(abs.(A), δ, ReachabilityAnalysis.BaseExp, Φ)
function _Φ₂_inv(A::AbstractMatrix, δ, alg, Φ=nothing)
    Aδ = A * δ
    if isnothing(Φ)
        Φ = _exp(Aδ, alg)
    end
    n = size(A, 1)
    N = eltype(A)
    In = Matrix(one(N) * I, n, n)
    B = Φ - In - Aδ
    Ainv = inv(Matrix(A))
    Ainvsqr = Ainv^2
    return Ainvsqr * B
end

# special case when A is a multiple of the identity
function _Φ₂_inv(A::IdentityMultiple, δ, alg, Φ=nothing)
    λ = A.M.λ
    @assert !iszero(λ) "the given identity multiple is not invertible"
    δλ = δ * λ
    α = (1 / λ)^2 * (exp(δλ) - 1 - δλ)
    return IdentityMultiple(α, size(A, 1))
end

function _Eplus(A::SparseMatrixCSC{N,D}, X0::LazySet{N}, δt; m=min(30, size(A, 1)),
                tol=1e-7) where {N,D}
    n = dim(X0)
    A2 = A * A # fast if A sparse
    V = symmetric_interval_hull(A2 * X0)
    v = V.radius
    Aabs = abs.(A)

    require(@__MODULE__, :ExponentialUtilities)
    Pv = _phiv(Aabs, v, 1, δt; m, tol)
    return Hyperrectangle(zeros(n), Pv)
end

# ---------------
# Absolute values
# ---------------

elementwise_abs(A::AbstractMatrix) = abs.(A)
function elementwise_abs(A::SparseMatrixCSC)
    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, abs.(nonzeros(A)))
end
elementwise_abs(A::IdentityMultiple) = IdentityMultiple(abs(A.M.λ), size(A, 1))

# =====================
# Optional dependencies
# =====================

using Requires: @require

function __init__()
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" eval(load_expokit_pade())
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" include("init_ExponentialUtilities.jl")
end

end  # module
