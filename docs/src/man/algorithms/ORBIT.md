```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Singleton propagation (ORBIT)

## Method

Consider the ODE

```math
x'(t) = Ax(t) + u, \qquad x_0 ∈ \mathbb{R}^n,\qquad (1)
```
with ``A ∈ \mathbb{R}^{n\times n}`` and ``u ∈ \mathbb{R}^n``. The analytic solution at ``δ > 0`` is known to be

```math
x(δ) = e^{Aδ}x_0 + ∫_0^δ e^{A(δ - s)} uds.
```
Let us introduce the matrices ``Φ(A, δ) := e^{Aδ}`` and ``Φ_1(A, δ) := ∫_0^{δ} e^{A(δ - s)} ds``. The computation of these matrices will be discussed later.

Since ``u`` is assumed constant,
```math
x(δ) = Φ x_0 + v,
```
where we have defined ``v := Φ_1(A, δ) u``. To get the exact solution at ``t = 2δ``, note that by time invariance of Eq. (1),

```math
x(2δ) = e^{Aδ} x(δ) + ∫_0^δ e^{A(δ - s)}u ds,
```
hence
```math
x(2δ) = Φ (Φ x_0 + v) + v = Φ^2 x_0 + Φ v + v.
```

The solution at subsequent multiples of the step-size ``δ`` is achieved by applying the preceding rule iteratlvely. For any ``k ≥ 1``, we obtain
```math
x(kδ) = Φ^k x_0 + \sum_{i=0}^{k-1} Φ^i v,\qquad k ≥ 1.
```

Before considering the computation of ``Φ`` and ``Φ_1``, note that the method admits a straightforward generalization for non-constant inputs. Indeed, if ``u : [0, T]\to \mathbb{R}^n`` is a piecewise-constant function, there always exists a sufficiently small step-size ``δ`` such that the following premise holds. Assume that ``\{u_1, u_2, …, u_N\} ⊆ \mathbb{R}^n`` is the range of values of the input, where ``u(t) = u_k`` for ``t ∈ [(k-1)δ, kδ)``, ``k = 1, …, N``, and let ``v_k := Φ_1(u_k, δ)`` for all ``k = 1,…, N``. Then it holds that ``x(δ) = Φ x_0 + v_{1}``, ``x(2δ) = Φ^2 x_0 + Φ v_1+ v_{2}``, and

```math
x(kδ) = Φ^k x_0 + \sum_{i=0}^{k-1} Φ^i v_{k-i},\qquad k = 1,…, N
```

The matrix ``Φ = e^{Aδ}`` can be evaluated in different ways, using the function
[`ReachabilityAnalysis.Exponentiation._exp`](@ref):

(1) `method = :base` uses Julia's built-in implementation (if `method = :base`),

(2) `method = :lazy` uses a lazy wrapper of the matrix exponential which is then evaluated using Krylov subspace methods.

Method (1) is the default method. Method (2) is particularly useful to work with
very large and sparse matrices (e.g. typically of order `n > 2000`). Evaluation
of ``Φ_1(u, δ)`` is available through the function
[`ReachabilityAnalysis.Exponentiation.Φ₁`](@ref). Two implementations are available:

(1) If the coefficients matrix ``A`` is invertible, then the integral is equivalent to computing ``A^{-1}(e^{Aδ} - I)``.

(2) In general, ``Φ_1(u, δ)`` can be obtained as a sub-block of a larger matrix. See [[FRE11]](@ref) for details.

Method (2) is the default method, although there are cases in which method (1) is more convenient.
For example, if we are only interested in singleton inputs and ``A`` is invertible,
it is possible to compute ``Φ_1(A, δ) u `` efficiently without actually inverting the matrix ``A`` in full.
