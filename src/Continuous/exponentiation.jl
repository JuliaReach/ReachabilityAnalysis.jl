using LinearAlgebra: checksquare
using LazySets.Arrays: isinvertible

# ==========================
# Exponentiation functions
# ==========================

# abstract supertype for all exponentiation methods
#abstract type AbstractExpMethod end

# matrix exponential using the scaling and squaring algorithim implemented in Julia Base
# :base
#struct BaseExp <: AbstractExpMethod end

# lazy wrapper for the matrix exponential, operations defined in LazySets.jl
# :lazy
#struct LazyExp <: AbstractExpMethod end

# matrix exponential for sparse matrices using Pade approximants, requires Expokit.jl
# :pade
#struct PadeExp <: AbstractExpMethod end

# general case: convert to Matrix
@inline _exp(A::AbstractMatrix, ::Val{:base}) = exp(Matrix(A))
@inline _exp(A::Matrix, ::Val{:base}) = exp(A)

# static arrays have their own exp method
@inline _exp(A::StaticArray, ::Val{:base}) = exp(A)

# exponential of an identity multiple, defined in MathematicalSystems.jl
@inline _exp(A::IdentityMultiple, ::Val{:base}) = exp(A)

# lazy wrapper of the matrix exponential
@inline _exp(A::AbstractMatrix, ::Val{:lazy}) = SparseMatrixExp(A)

# pade approximants (requires Expokit.jl)
@inline _exp(A::SparseMatrixCSC, ::Val{:pade}) = padm(A)

"""
    _exp(A::AbstractMatrix, δ, method)

Compute the matrix exponential ``e^{Aδ}``.

### Input

- `A`       -- matrix
- `δ`       -- step size
- `method`  -- the method used to take the matrix exponential of `A`;
               possible options are:

    - `:base` -- use Higham's scaling and squaring method implemented in Julia
                 standard library; see `?exp` for details
    - `:lazy` -- return a lazy wrapper type around the matrix exponential using
                 the implementation `LazySets.SparseMatrixExp`
    - `:pade` -- apply Pade approximant method to compute the matrix exponential
                 of a sparse matrix (requires `Expokit`)

### Output

A matrix or a lazy wrapper of the matrix exponential, depending on `method`.

### Notes

If the algorithm `"lazy"` is used, evaluations of the action of the matrix
exponential are done with the `expmv` implementation from `Expokit`
(but see `LazySets#1312` for the planned generalization to other backends).
"""
function _exp(A::AbstractMatrix, δ, method)
    n = checksquare(A)
    return _exp(A * δ, method)
end

function _exp(A::AbstractMatrix, δ, method::Val{:pade})
    n = checksquare(A)
    @requires Expokit
    return _exp(A * δ, method)
end

# default
_exp(A::AbstractMatrix, δ) = _exp(A, δ, Val(:base))

@inline function _Aδ_3n(A::AbstractMatrix{N}, δ, n) where {N}
    return [Matrix(A*δ)     Matrix(δ*I, n, n)  zeros(n, n)    ;
            zeros(n, 2*n          )  Matrix(δ*I, n, n)        ;
            zeros(n, 3*n          )                           ]::Matrix{N}
end

@inline function _Aδ_3n(A::SparseMatrixCSC{N, M}, δ, n) where {N, M}
    return [sparse(A*δ)     sparse(δ*I, n, n)  spzeros(n, n) ;
            spzeros(n, 2*n          )  sparse(δ*I, n, n)     ;
            spzeros(n, 3*n          )                        ]::SparseMatrixCSC{N, M}
end

"""
    Φ₁(A, δ, [method]=:base, [isinv]=false)

Evaluate the series

```math
Φ₁(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+1}}{(i+1)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`      -- coefficients matrix
- `δ`      -- step size
- `method` -- (optional, default: `:base`) the method used to take the matrix
              exponential of the coefficient matrix; see the documentation of
              [`_exp`](@ref) for available options
- `isinv`  -- (optional, default: `false`) if `true`, assume that the coefficients
              matrix is invertible

### Output

A matrix.

### Algorithm

If ``A`` is invertible, ``Φ₁`` can be computed as

```math
Φ₁(A, δ) = A^{-1}(e^{δA} - I_n),
```
where ``I_n`` is the identity matrix of order ``n``.

In the general case, implemented in this function, it can be
computed as submatrices of the block matrix

```math
P = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\\\
0 && 0     && δI_n \\\\
0 && 0     && 0
\\end{pmatrix}.
```
It can be shown that `Φ₁(A, δ) = P[1:n, (n+1):2*n]`.
We refer to [[FRE11]](@ref) for details.
"""
function Φ₁(A::AbstractMatrix, δ, method=:base, isinv=false)
    return _Φ₁(A, δ, Val(method), Val(isinv))
end

# dispatch for discretization methods
__Φ₁(A, δ, method, isinv::Val{:false}, Φ, U) = _Φ₁(A, δ, method, isinv)
__Φ₁(A, δ, method, isinv::Val{:true}, Φ, U) = _Φ₁(A, δ, method, isinv, Φ)
__Φ₁(A, δ, method, isinv::Val{:true}, Φ, U::Singleton) = _Φ₁(A, δ, method, isinv, U, Φ)

function _Φ₁(A::AbstractMatrix, δ, method, isinv::Val{:false})
    n = checksquare(A)
    B = _Aδ_3n(A, δ, n)
    P = _exp(B, method)
    return _P₁_blk(P, n, method)
end

@inline _P₁_blk(P, n::Int, ::Val{:base}) = P[1:n, (n+1):2*n]
@inline _P₁_blk(P, n::Int, ::Val{:lazy}) = sparse(get_columns(P, (n+1):2*n)[1:n, :])
@inline _P₁_blk(P, n::Int, ::Val{:pade}) = P[1:n, (n+1):2*n]

# compute the matrix Φ₁(A, δ) = A^{-1}(exp(Aδ) - I), assuming that A is invertible
function _Φ₁(A::AbstractMatrix{N}, δ, method, isinv::Val{:true}, Φ=nothing) where {N}
    if isnothing(Φ)
        Φ = _exp(A, δ, method)
    end
    n = size(A, 1)
    In = Matrix(one(N)*I, n, n)
    Ainv = inv(A)
    return Ainv * (Φ - In)
end

# compute the vector Φ₁(A, δ)U = A^{-1}(exp(Aδ) - I) U, where U is a vector,
# assuming that A is invertible
function _Φ₁(A::AbstractMatrix, δ, method, isinv::Val{:true}, U::AbstractVector, Φ=nothing)
    if isnothing(Φ)
        Φ = _exp(A, δ, method)
    end
    x = Φ * U - U
    return A \ x
end

function load_Φ₁_krylov()
return quote

function _Φ₁(A::AbstractMatrix, δ, method::Val{:lazy}, isinv::Val{:true}, U::AbstractVector; m=30, tol=1e-10)
    w = expv(1.0 * δ, A, U, m=m, tol=tol)
    x = w - U
    return A \ x
end

end end  # quote / load_Φ₁_krylov()

"""
    Φ₂(A::AbstractMatrix, δ, [method]=:base, [isinv]=false)

Evaluate the series

```math
Φ₂(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+2}}{(i+2)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`      -- coefficient matrix
- `δ`      -- step size
- `method` -- (optional, default: `:base`) the method used to take the matrix
              exponential of the coefficient matrix; see the documentation of
              [`_exp`](@ref) for available options
- `isinv`  -- (optional, default: `false`) if `true`, assume that the coefficients
              matrix is invertible

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
P = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\\\
0 && 0 && δI_n \\\\
0 && 0 && 0
\\end{pmatrix}.
```
It can be shown that `Φ₂(A, δ) = P[1:n, (2*n+1):3*n]`.
We refer to [[FRE11]](@ref) for details.
"""
function Φ₂(A::AbstractMatrix, δ, method=:base, isinv=false)
    return _Φ₂(A, δ, Val(method), Val(isinv))
end

# dispatch for discretization methods
__Φ₂(A, δ, method, isinv::Val{:false}) = _Φ₂(A, δ, method, isinv)
__Φ₂(A, δ, method, isinv::Val{:false}, Φ) = _Φ₂(A, δ, method, isinv)

__Φ₂(A, δ, method, isinv::Val{:true}, Φ) = _Φ₂(A, δ, method, isinv, Φ)
__Φ₂(A, δ, method, isinv::Val{:true}) = _Φ₂(A, δ, method, isinv)

function _Φ₂(A::AbstractMatrix, δ, method, isinv::Val{:false})
    n = checksquare(A)
    B = _Aδ_3n(A, δ, n)
    P = _exp(B, method)
    return _P₂_blk(P, n, method)
end

@inline _P₂_blk(P, n, ::Val{:base}) = P[1:n, (2*n+1):3*n]
@inline _P₂_blk(P, n, ::Val{:lazy}) = sparse(get_columns(P, (2*n+1):3*n)[1:n, :])
@inline _P₂_blk(P, n, ::Val{:pade}) = P[1:n, (2*n+1):3*n]

# compute the matrix Φ₂ = A^{-2} (exp(A*δ) - I - A*δ) assuing A is invertible
function _Φ₂(A::AbstractMatrix{N}, δ, method, isinv::Val{:true}, Φ=nothing) where {N}
    # @assert isinvertible(A) "the given matrix should be invertible"
    if isnothing(Φ)
        Φ = _exp(A, δ, method)
    end
    n = size(A, 1)
    In = Matrix(one(N)*I, n, n)
    B = Φ - In - Aδ
    Ainv = inv(A)
    Ainvsqr = Ainv^2
    return Ainvsqr * B
end

function _Φ₂(A::IdentityMultiple, δ, method, isinv::Val{:true}, Φ=nothing)
    λ = A.M.λ
    @assert !iszero(λ) "the given identity multiple is not invertible"
    δλ = δ * λ
    α = (1/λ)^2 * (exp(δλ) - 1 - δλ)
    return IdentityMultiple(α, size(A, 1))
end

# ================
# Absolute values
# ================

@inline _elementwise_abs(A::AbstractMatrix) = abs.(A)
@inline _elementwise_abs(A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, abs.(nonzeros(A)))
@inline _elementwise_abs(A::IdentityMultiple) = IdentityMultiple(exp(A.M.λ), size(A, 1))

# ====================================
# Exponentiation of interval matrices
# ====================================

# TODO: outsource to IntervalMatrices.jl

# compute Iδ + 1/2 * δ^2 * A + 1/6 * δ^3 * A² + ... + 1/(η+1)! * δ^(η+1) * A^η + E(δ) * δ
function _Cδ(A, δ, order)
    n = size(A, 1)
    A² = A * A
    if isa(A, IntervalMatrix)
        Iδ = IntervalMatrix(Diagonal(fill(IA.Interval(δ), n)))
    else
        Iδ = Matrix(δ * I, n, n)
    end

    IδW = Iδ + 1/2 * δ^2 * A + 1/6 * δ^3 * A²
    M = IδW

    if order > 2
        # i = 2
        αᵢ₊₁ = 6 # factorial of (i+1)
        Aⁱ = A²
        δⁱ⁺¹ = δ^3
        @inbounds for i in 3:order
            αᵢ₊₁ *= i+1
            δⁱ⁺¹ *= δ
            Aⁱ *= A
            M += (δⁱ⁺¹/αᵢ₊₁) * Aⁱ
        end
    end
    E = IntervalMatrices._exp_remainder(A, δ, order; n=n)
    return M + E * δ
end
