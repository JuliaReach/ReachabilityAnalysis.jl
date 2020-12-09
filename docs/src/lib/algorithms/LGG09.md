```@docs
LGG09
```

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

Let us write the recurrence

```math
X_{k+1} = \Phi X_k \oplus V_k,\qquad 0 \leq k \leq N
```
in the unrapped form,

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

## Parallel formulation

The support functions of the sequence $\{X_k\}_k$ along different directions can be
computed in parallel. Indeed, if $d_1$ and $d_2$ are two given template directions, two different processes
can independently compute $\rho(d_1, X_k)$ and $\rho(d_2, X_k)$ for all $k = 0, 1, \ldots, N$
using the methods described above. Once both computations have finished, we can store
the resulting support functions in the same array. Use the flag `threaded=true` to
use this method.

Implementation-wise the function `_reach_homog_dir_LGG09!` spawns differen threads
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
In more generality, it holds that $(\Phi^T)^k d  = \lambda^k d$.
Applying this relation to the support function recurrence described above, we get
for the general inhomogeneous and possibly time-varying inputs case:

```math
\rho(d, X_k) = \rho(\lambda^k d, X_0) + \sum_{i=0}^{k-1} \rho(\lambda^{k-i-1} d, V_i).
```
To further simplify this formula, we analyze different cases of $\lambda$.
If $\lambda = 0$, then $\rho(d, X_k) = \rho(d, V_k)$ for all $k \geq 0$, so we focus
on either $\lambda$ being positive or negative. To further simplify the computation
of $\rho(d, X_k)$, we can use the property $\rho(\lambda d, X) = \lambda \rho(d, X)$
if $\lambda \geq 0$. We now consider the cases $\lambda > 0$ and $\lambda < 0$.

**Case $\lambda > 0$.** Then $\lambda^k > 0$ for all $k \geq 1$,
and

```math
\rho(d, X_k) = \lambda^k \rho(d, X_0) +  \sum_{i=0}^{k-1} \lambda^{k-i-1} \rho(d, V_i).
```
In the particular case that $V$ is constant, we can extract $\rho(d, V)$ from the sum
and are left with evaluating the support function only at $\rho(d, X_0)$  and $\rho(d, V)$
to construct the full sequence $\{\rho(d, X_k)\}_{k}$.

**Case $\lambda < 0$.** Then $\lambda^k = (-1)^k (-\lambda)^k$, so
$\lambda^k$ is positive if $k$ is even and it is negative otherwise. Note that

```math
\rho(d, X_k) = (-\lambda)^k \rho((-1)^k d, X_0) + \sum_{i=0}^{k-1} (-\lambda)^{k-i-1} \rho((-1)^{k-i-1} d, V_i).
```
In this case, if $V$ is constant we have to compute four support functions, namely
$\rho(\pm d, X_0)$ and $\rho(\pm d, V)$.

## Implementation details

The computed support function values can accessed directly through the field
`sf::SN` of each template reach-set. Here `sf` is an array view of `ρℓ::Matrix{N}(undef, length(dirs), NSTEPS)`;
note that each row corresponds to a template direction and each column corresponds to a fixed
iteration index $k \geq 0$.
