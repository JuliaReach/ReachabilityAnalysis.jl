```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Support-function based method (LGG09)

## Method

In the following subsections we outline the method of [[LGG09]](@ref) to solve
linear set-based recurrences using support functions, first the homogeneous case
and then the inhomogeneous case without wrapping effect.
We also discuss the special case of real eigenvalues.

## Homogeneous case

Consider the set-based recurrence

```math
X_{k+1} = \Phi X_k,\qquad 0 \leq k \leq N
```
where $\Phi \in \mathbb{R}^{n\times n}$ and $X_0 \subseteq \mathbb{R}^n$ are given.
By unrwapping the recurrence, $X_k = \Phi^k X_0$ for all $k \geq 0$. Let $d \in \mathbb{R}^n$ be a given
*template direction*. Using the property of support functions $\rho(d, A X) = \rho(A^T d, X)$
for any matrix $A$ and set $X$, we have that

```math
ρ(d, X_k) = \rho(d, \Phi^k X_0) = \rho((\Phi^T)^k d, X_0).
```
In this way we are able to reason with the sequence $\{X_0, X_1, X_2, \ldots, X_N\}$
by evaluating the support function of the initial set $X_0$ along the directions
$\{d, \Phi^T d, (\Phi^T)^2 d, \ldots, (\Phi^T)^N d\}$.


## Inhomogeneous case

The inhomogeneous case generalizes the previous case by taking, at each step,
the Minkowski sum with an element from the sequence $\{V_0, V_1, V_2, \ldots, V_N\}$:

```math
X_{k+1} = \Phi X_k \oplus V_k,\qquad 0 \leq k \leq N.
```
Let us write such recurrence in the unrapped form,

```math
\begin{aligned}
\quad X_1 &= \Phi X_0 \oplus V_0 \\[1mm]
\quad X_2 &= \Phi X_1 \oplus V_1 = \Phi^2 X_0 \oplus \Phi V_0 \oplus V_1 \\[1mm]
\quad X_3 &= \Phi X_2 \oplus V_2 = \Phi^3 X_0 \oplus \Phi^2 V_0 \oplus \Phi V_1  \oplus V_2 \\[1mm]
\quad &\vdots \\[1mm]
\quad X_k &= \Phi^k X_0 \oplus \left( \bigoplus_{i=0}^{k-1} \Phi^{k-i-1} V_i \right)
\end{aligned}
```
where the big Minkowski sum is just an abbreviation for
$\Phi^{k-1} V_0 \oplus \Phi^{k-2} V_1 \oplus \Phi^{k-3} V_2 \oplus \ldots \oplus \Phi V_{k-2} \oplus V_{k-1}$.

Let $d \in \mathbb{R}^n$ be a given template direction. Using the additive property of
support functions, $\rho(d, X \oplus Y) = \rho(d, X) + \rho(d, Y)$ for any sets $X$ and $Y$,
we have that

```math
\begin{aligned}
\quad \rho(d, X_1) &= \rho(\Phi^T d, X_0) + \rho(d, V_0) \\[1mm]
\quad \rho(d, X_2) &= \rho((\Phi^T)^2 d, X_0) + \rho(\Phi^T d, V_0) + \rho(d, V_1) \\[1mm]
\quad \rho(d, X_3) &= \rho((\Phi^T)^3 d, X_0) + \rho((\Phi^T)^2 d, V_0) + \rho(\Phi^T d, V_1) + \rho(d, V_2) \\[1mm]
\quad &\vdots \\[1mm]
\quad \rho(d, X_k) &= \rho((\Phi^T)^k d, X_0) + \sum_{i=0}^{k-1} \rho( (\Phi^T)^{k-i-1} d,  V_i).
\end{aligned}
```
In a similar fashion to the homogeneous case, the method allows to efficiently reason
about the the sequence $\{X_0, X_1, X_2, \ldots, X_N\}$ by evaluating the support
function of the initial set $X_0$ and the input sets $\{V_k\}_k$ along the directions
$\{d, \Phi^T d, (\Phi^T)^2 d, \ldots, (\Phi^T)^N d\}$. Implementation-wise, we update
two sequences, one that accounts for the homogeneous term, and another
sequence that accounts for the effect of the accumulated inputs.

## Implementation details

The reach-set representation used is a [`TemplateReachSet`](@ref), which stores the
directions used (vector of vectors) and the support function evaluated at each direction
(matrix, see below). The set representation, `set(R::TemplateReachSet)`, is either a polyhedron in constraint form
(`HPolyhedron`), or a polytope (`HPolytope`) if the directions are bounding, i.e.
the template directions define a bounded set.

The computed support function values can accessed directly through the field
`sf::SN` of each template reach-set. Here `sf` is an array view of type `::Matrix{N}(undef, length(dirs), NSTEPS)`:
each row corresponds to one of the template directions and each column corresponds to a fixed iteration index $k \geq 0$.

If you use directions from the canonical basis of $\mathbb{R}^n$, it is recommended to define `LazySets.Arrays.SingleEntryVector`
or "one-hot" arrays as they are commonly called, because there are methods that dispatch on such type of arrays efficiently.

## Parallel formulation

The support functions of the sequence $\{X_k\}_k$ along different directions can be
computed in parallel. Indeed, if $d_1$ and $d_2$ are two given template directions, two different processes
can independently compute $\rho(d_1, X_k)$ and $\rho(d_2, X_k)$ for all $k = 0, 1, \ldots, N$
using the methods described above. Once both computations have finished, we can store
the resulting support functions in the same array. Use the flag `threaded=true` to
use this method.

Implementation-wise the function `_reach_homog_dir_LGG09!` spawns different threads
which populate the matrix `ρℓ::Matrix{N}(undef, length(dirs), NSTEPS)` with the computed
values. Hence each thread computes a subset of distinct rows of `ρℓ`.

## Real eigenvalues

If the spectrum of the state matrix only has real eigenvalues, the sequence of
support functions can be computed efficiently if we work with a template
consisting of eigenvectors of $\Phi^T$. This idea is described in [[LGG09]](@ref)
and we recall it here for convenience.

The method stems from the fact that if $(\lambda, d)$ is an eigenvalue-eigenvector
pair of the matrix $\Phi^T$, with $\lambda \in \mathbb{R}$, then
$\Phi^T d = \lambda d$, and if we apply $\Phi^T$ on both sides of this identity, we get
$(\Phi^T)^2 d = \Phi^T (\Phi^T d) = \Phi^T(\lambda d) = \lambda^2 d$.
In more generality, it holds that $(\Phi^T)^k d  = \lambda^k d$ for all $k \ge 1$.
Applying this relation to the support function recurrence described above, we get
for the general inhomogeneous and possibly time-varying inputs case:

```math
\rho(d, X_k) = \rho(\lambda^k d, X_0) + \sum_{i=0}^{k-1} \rho(\lambda^{k-i-1} d, V_i).
```
To further simplify this formula, we analyze different cases of $\lambda$.
If $\lambda = 0$, then $\rho(d, X_k) = \rho(d, V_k)$ for all $k \geq 1$, so we focus
on either $\lambda$ being positive or negative. To further simplify the computation
of $\rho(d, X_k)$, we can use the property $\rho(\lambda d, X) = \lambda \rho(d, X)$
if $\lambda \geq 0$. We now consider the cases $\lambda > 0$ and $\lambda < 0$.

**Case $\lambda > 0$.** Then $\lambda^k > 0$ for all $k \geq 1$,
and

```math
\rho(d, X_k) = \lambda^k \rho(d, X_0) +  \sum_{i=0}^{k-1} \lambda^{k-i-1} \rho(d, V_i).
```
We are left with evaluating the support function only at $\rho(d, X_0)$  and $\rho(d, V_i)$
to construct the full sequence $\{\rho(d, X_k)\}_{k}$. Moreover, if the $V_i$'s are constant
we can extract them from the right-hand side sum and use that
```math
\sum_{i=0}^{k-1} \lambda^{k-i-1} = 1 + \lambda + \ldots + \lambda^{k-1} = \dfrac{1 - \lambda^k}{1 - \lambda}.
```

**Case $\lambda < 0$.** Since $\lambda^k = (-1)^k (-\lambda)^k$ and $\lambda < 0$, then
$\lambda^k$ is positive if $k$ is even, otherwise it is negative. So we can write:

```math
\rho(d, X_k) = (-\lambda)^k \rho((-1)^k d, X_0) + \sum_{i=0}^{k-1} (-\lambda)^{k-i-1} \rho((-1)^{k-i-1} d, V_i).
```
The main difference between this case and the previous one is that now we have to evaluate
support functions $\rho(\pm d, X_0)$ and $\rho(\pm d, V_i)$. Again, simplification takes place
if the $V_i$'s are constant and such special case is considered in the implementation.
