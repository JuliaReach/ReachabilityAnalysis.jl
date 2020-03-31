using LinearAlgebra: checksquare
using LazySets.Arrays: isinvertible

# ==========================
# Exponentiation functions
# ==========================

# general case: convert to Matrix
@inline _exp(A::AbstractMatrix) = exp(Matrix(A))

# static arrays have their own exp method
@inline _exp(A::StaticArray) = exp(A)

# lazy wrapper (defined in LazySets)
@inline _exp_lazy(A::AbstractMatrix) = SparseMatrixExp(A)

# pade approximants (requires Expokit.jl)
@inline _exp_pade(A::SparseMatrixCSC) = padm(A)

# exponential of an identity multiple, defined in MathematicalSystems
@inline _exp(A::IdentityMultiple) = exp(A)

"""
    _exp(A::AbstractMatrix, δ::Float64, method::Symbol)

Compute the matrix exponential ``e^{Aδ}``.

### Input

- `A`       -- matrix
- `δ`       -- step size
- `method`  -- symbol with the method used to take the matrix exponential of `A`;
               possible options are:

    - `:base` -- use the scaling and squaring method implemented in Julia standard
                 library; see `?exp` for details
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
function _exp(A::AbstractMatrix, δ::Float64, method::Symbol)
    n = checksquare(A)
    if method == :base
        return _exp(A * δ) # TODO use dots ? (requires MathematicalSystems#189 for IdentityMultiple)

    elseif method == :lazy
        return _exp_lazy(A * δ)

    elseif method == :pade
        @requires Expokit
        return _exp_pade(A * δ)

    else
       throw(ArgumentError("the exponentiation method $method is unknown"))
    end
end

@inline function _Aδ_3n(A::AbstractMatrix, δ::Float64, n::Int)
    return [A*δ     sparse(δ*I, n, n)  spzeros(n, n)    ;
            spzeros(n, 2*n          )  sparse(δ*I, n, n);
            spzeros(n, 3*n          )                   ]
end

"""
    Φ₁(A, δ, method)

Compute the series

```math
Φ₁(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+1}}{(i+1)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
- `method`  -- (optional, default: `"base"`) the method used to take the matrix
                    exponential of the coefficient matrix; see the documentation of
                    [`exp_Aδ`](@ref) for available options

### Output

A matrix.

### Algorithm

We use the method from [1]. If ``A`` is invertible, ``Φ₁`` can be computed as

```math
Φ₁(A, δ) = A^{-1}(e^{δA} - I_n),
```
where ``I_n`` is the identity matrix of order ``n``.

In the general case, implemented in this function, it can be computed as
submatrices of the block matrix

```math
P = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\\\
0 && 0     && δI_n \\\\
0 && 0     && 0
\\end{pmatrix}.
```
It can be shown that `Φ₁(A, δ) = P[1:n, (n+1):2*n]`.

[1] Frehse, Goran, et al. "SpaceEx: Scalable verification of hybrid systems."
International Conference on Computer Aided Verification. Springer, Berlin,
Heidelberg, 2011.
"""
function Φ₁(A::AbstractMatrix, δ::Float64, method::String)
    n = checksquare(A)
    B = _Aδ_3n(A, δ, n)

    if method == :base
        P = _exp(B)
        return P[1:n, (n+1):2*n]

    elseif method == :lazy
        P = _exp_lazy(B)
        return sparse(get_columns(P, (n+1):2*n)[1:n, :])

    elseif method == :pade
        P = _exp_pade(B)
        return P[1:n, (n+1):2*n]

    else
        throw(ArgumentError("the exponentiation method $method is unknown"))
    end
end

"""
    Φ₂(A, δ, method)

Compute the series

```math
Φ₂(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+2}}{(i+2)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
- `method`       -- the method used to take the matrix
                    exponential of the coefficient matrix; see the documentation of
                    [`exp_Aδ`](@ref) for available options

### Output

A matrix.

### Algorithm

We use the method from [1]. If ``A`` is invertible, ``Φ₂`` can be computed as

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

[1] Frehse, Goran, et al. "SpaceEx: Scalable verification of hybrid systems."
International Conference on Computer Aided Verification. Springer, Berlin,
Heidelberg, 2011.
"""
function Φ₂(A::AbstractMatrix, δ::Float64, method::Symbol)
    if method == :inverse
        return _Φ₂_inverse(A, δ)
    end

    n = checksquare(A)
    B = _Aδ_3n(A, δ, n)

    if method == :base
        P = _exp(B)
        return P[1:n, (2*n+1):3*n]

    elseif method == :lazy
        P = _exp_lazy(B)
        return sparse(get_columns(P, (2*n+1):3*n)[1:n, :])

    elseif method == :pade
        P = _exp_pade(B)
        return P[1:n, (2*n+1):3*n]

    else
       throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end
end

@inline function _Φ₂_inverse(A::IdentityMultiple, δ::Float64)
    λ = A.M.λ
    @assert !iszero(λ) "the given identity multiple is not invertible"
    δλ = δ * λ
    α = (1/λ)^2 * (exp(δλ) - 1 - δλ)
    return IdentityMultiple(α, size(A, 1))
end

@inline function _Φ₂_inverse(A::AbstractMatrix, δ::Float64)
    @assert isinvertible(A) "the given matrix should be invertible"
    Ainv = inv(A)
    n = size(A, 1)
    In = Matrix(one(N)*I, n, n)
    return Ainv^2 * (exp(δ*A) - In - δ * A)
end

# ================
# Absolute values
# ================

@inline _elementwise_abs(A::AbstractMatrix) = abs.(A)
@inline _elementwise_abs(A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, abs.(nonzeros(A)))
@inline _elementwise_abs(A::IdentityMultiple) = IdentityMultiple(exp(A.M.λ), size(A, 1))
