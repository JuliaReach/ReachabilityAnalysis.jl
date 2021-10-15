```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Singleton propagation (ORBIT)

## Method

Consider the ODE

```math
x'(t) = Ax(t) + u, \qquad x_0 \in \mathbb{R}^n,\qquad (1)
```
with $A \in \mathbb{R}^{n\times n}$ and $u \in \mathbb{R}^n$. The analytic solution at $\delta > 0$ is known to be

```math
x(\delta) = e^{A\delta}x_0 + \int_0^\delta e^{A(\delta - s)} uds.
```
Let us introduce the matrices $\Phi(A, \delta) := e^{A\delta}$ and $\Phi_1(A, \delta) := \int_0^{\delta} e^{A(\delta - s)} ds$. The computation of these matrices will be discussed later.

Since $u$ is assumed constant,
```math
x(\delta) = \Phi x_0 + v,
```
where we have defined $v := \Phi_1(A, \delta) u$. To get the exact solution at $t = 2\delta$, note that by time invariance of Eq. (1),

```math
x(2\delta) = e^{A\delta} x(\delta) + \int_0^\delta e^{A(\delta - s)}u ds,
```
hence
```math
x(2\delta) = \Phi (\Phi x_0 + v) + v = \Phi^2 x_0 + \Phi v + v.
```

The solution at subsequent multiples of the step-size $\delta$ is achieved by applying the preceding rule iteratlvely. For any $k \geq 1$, we obtain
```math
x(k\delta) = \Phi^k x_0 + \sum_{i=0}^{k-1} \Phi^i v,\qquad k \geq 1.
```

Before considering the computation of $\Phi$ and $\Phi_1$, note that the method admits a straightforward generalization for non-constant inputs. Indeed, if $u : [0, T]\to \mathbb{R}^n$ is a piecewise-constant function, there always exists a sufficiently small step-size $\delta$ such that the following premise holds. Assume that $\{u_1, u_2, \ldots, u_N\} \subseteq \mathbb{R}^n$ is the range of values of the input, where $u(t) = u_k$ for $t \in [(k-1)\delta, k\delta)$, $k = 1, \ldots, N$, and let $v_k := \Phi_1(u_k, \delta)$ for all $k = 1,\ldots, N$. Then it holds that $x(\delta) = \Phi x_0 + v_{1}$, $x(2\delta) = \Phi^2 x_0 + \Phi v_1+ v_{2}$, and

```math
x(k\delta) = \Phi^k x_0 + \sum_{i=0}^{k-1} \Phi^i v_{k-i},\qquad k = 1,\ldots, N
```

The matrix $\Phi = e^{A\delta}$ can be evaluated in different ways, using the function [`_exp`](@ref):

(1) `method=:base` uses Julia's built-in implementation (if `method = :base`),

(2) `method = :lazy` uses a lazy wrapper of the matrix exponential which is then evaluted using Krylov subspace methods.

Method (1) is the default method. Method (2) is particularly useful to work with very large and sparse matrices (e.g. typically of order `n > 2000`). Evaluation of $\Phi_1(u, \delta)$ is available through the function [`Φ₁`](@ref). Two implementations are available:

(1) If the coefficients matrix $A$ is invertible, then the integral is equivalent to computing $A^{-1}(e^{A\delta} - I)$.

(2) In general, $\Phi_1(u, \delta)$ can be obtained as a sub-block of a larger matrix. See [[FRE11]](@ref) for details.

Method (2) is the default method, although there are cases in which method (1) is more convenient.
For example, if we are only interested in singleton inputs and $A$ is invertible,
it is possible to compute $\Phi_1(A, \delta) u $ efficiently without actually inverting the matrix $A$ in full.
