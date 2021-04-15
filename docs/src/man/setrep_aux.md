### FULL SETREP STUFF


```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

# Set representations

Subets of $\mathbb{R}^n$ can be represented in many different ways, and each
reprsentation has its advantages and disadvantages: there is no one-size-fits-all
set represententation in reachability applications. Depending on the type of
operation that we want to apply, one set representation may be more convenient
than another one.

Two important characteristics are:

- **Closure:** A set representation is said to be *closed* under a given operation
  if the result of the operation is again a set in the same class.

- **Cost:** To measure the cost of making a given operation on a set, we consider
  the total number of binary operations, denoted with $\mathrm{Op}(\cdot)$. (....)

In the rest of this section we define the most common set representations used in
reachability analysis and recall some fundamental properties. Moreover, we show
examples using the [LazySets.jl](https://github.com/JuliaReach/LazySets.jl).


## Introduction

In simple terms, reachability analysis is concerned with studying the sets of states
that a system can reach, starting from a set of initial states and under the influence
of a set of input trajectories and parameter values.

To fix ideas, consider a purely continuous dynamical system over some state
space $\mathcal{X} \subset \mathbb{R}^n$ defined via a differential equation of the form:

```math
    x'(t) = f(x, u),
```
where $u(t) \in \mathcal{U}(t)$ ranges over some specified set of admissible input signals.
Given a set of initial states $\mathcal{X}_0 \subseteq \mathcal{X}$, the problem we are
interested in is to compute *all* the states visited by trajectories of the system starting
from any $x_0 \in \mathcal{X}_0$ and for any admissible input function $u(t)$. In other
sections of this manual we also consider generalizations of this problem such as uncertain
parameters, or dynamical systems which are hybrid, i.e. mixing continuous dynamics and discrete transitions.

In order to compute with sets, we need to find convenient ways to represent and operate with them.

In the next section we consider some common set representations and show how to work with `LazySets.jl`.

## Reach-sets

More technically, we define the reachable set at a given time point
$\delta \in \mathbb{R}$, also known as the *reach-set* for the ODE
```math
x'(t) = f(x(t), u(t), p(t), t),
```
as given by
```math
\mathcal{R}(δ) := \left\{ x(δ) = \int_0^δ f(x(t), u(t), p(t), t) dt, x(0) ∈ X_0, u(t) ∈ \mathcal{U}, p(t) ∈ \mathcal{P} \right\}.
```
Here $X_0$ denotes the set of initial states, $\mathcal{U}$ denotes the input set,
and $\mathcal{P}$ denotes the parameter values. For practical problems, the set
$\mathcal{R}(δ)$ cannot be obtained exactly, and reachability methods aim at
computing suitable over-approximations (or under-approximations) of it.

## Flowpipes



### Polyhedra

A hyperplane is the set $\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx = b\}$, where $a ∈ \mathbb{R}^n$ is the normal vector and $b ∈ \mathbb{R}$ is the displacement.

In `LazySets.jl`, the type `hyperplane` can be used to define a hyperplane. For example, . . . .

EJEMPLO <<<<<

A halfspace is the set $\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx ≤ b\}$, where $a ∈ \mathbb{R}^n$ is the normal vector and $b ∈ \mathbb{R}$ is the displacement. Please note that a half-space defines the region on one side of the hyperplane ``a^Tx = b``.

In `LazySets.jl`, the type `HalfSpace` can be used to define a half-space. For example, . . . .

EJEMPLO <<<<<

When we consider the intersection of a finite subset of half-spaces, we get a polyhedron. A polyhedron is thus the set defined as $\mathcal{P} ⊆ \mathbb{R}^n$,
```math
    \mathcal{P} = \{x \in \mathbb{R}^n | \bigcap_{i=1}^m a^T_i x \leq b_i \},
```
where $a_i \in \mathbb{R}^n$ and $b_i \in \mathbb{R}$. A polytope is a bounded polyhedron.

In `LazySets.jl`, the types `HPolyhedron` and `HPolytope` representent polyhedron and polytopes respectively. For example, . . . .



### Support functions

The support function of a compact set $\mathcal{X} \subseteq \mathbb{R}^n$ attributes to a direction $\ell \in \mathbb{R}^n$ the real number
```math
    \rho_{\mathcal{X}} (\ell) = \max\{ \ell^T x | x \in \mathcal{X} \}.
```
It is important to note that for a given direction $\ell$, the support function defines the position of a halfspace
```math
    \mathcal{H}_{\ell} = \{x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell)\}.
```

```@example
using LazySets, Plots

E = Ellipsoid([1.0, 1.0], [2.7 0.28; 0.28 1.04])

H1 = HalfSpace([1.0, 2.0], ρ([1.0, 2.0], E))
H2 = HalfSpace([-1.0, -2.0], ρ([-1.0, -2.0], E))
H3 = HalfSpace([-1.0, 2.0], ρ([-1.0, 2.0], E))
H4 = HalfSpace([1.0, -2.0], ρ([1.0, -2.0], E))


P = HPolytope([H1, H2, H3, H4])

plot(E)
plot!(P)

plot!(Singleton(σ([1.0, 2.0], E)), lab="d1")
plot!(Singleton(σ([-1.0, -2.0], E)), lab="-d1")
plot!(Singleton(σ([-1.0, 2.0], E)), lab="d2")
plot!(Singleton(σ([1.0, -2.0], E)), lab="-d2")
```

The following formulas hold:

- $\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell) =

- $\rho_{M\mathcal{X}}(\ell) =  for any
Here $M\mathcal{X}$ represents the linear map for $M \in \mathbb{R}^{m \times n}$.
- $\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell) =

The following table summarizes the number of
| Operation                                     | Simplification rule | Cost  |
|-----------------------------------------------|-------|------|
| $\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell)$ |$\rho_{\mathcal{X}}(\ell) + \rho_{\mathcal{Y}}(\ell)$  |  $1$     |
| $\rho_{M\mathcal{X}}(\ell)$                   |$rho_{\mathcal{X}}(M^T\ell)$|   $2mn$   |
| $\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell)$   |$\max\{\rho_{\mathcal{X}(\ell)}, \rho_{\mathcal{Y}(\ell)}\}$|  $1$    |

which touches and contains $\mathcal{X}$ . If $\ell$ is of unit length, then
$\rho_{\mathcal{X}}(\ell)$ is the signed distance of $\mathcal{H}_{\ell}$ to the origin.
Evaluating the support function for a set of directions $L ⊆ \mathbb{R}^n$ provides an overapproximation
```math
    \lceil \mathcal{X} \rceil _L = \bigcap_{\ell \in L} \{ x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell) \}
```

i.e., $\mathcal{X} ⊆ \lceil \mathcal{X} \rceil _L$. If $L = \mathbb{R}^n$, then $\mathcal{X} = \lceil \mathcal{X} \rceil _L$, so the support function represents any convex set $\mathcal{X}$ exactly. If $L$ is a finite set of directions $L = {\ell_1, . . . , \ell_m}$, then $\lceil \mathcal{X} \rceil _L$ is a polyhedron.

### Hyperrectangular sets

A special class of polyhedra are (...)


### Zonotopic sets

Zonotopes are a sub-class of polytopes defined as the image of a unit cube under
an affine transformation. An equivalent characterization of zonotopes is the
generator representation. Here, $Z ⊆ \mathbb{R}^n$ is defined by a center $c ∈ \mathbb{R}^n$ and a finite number of generators $g_1, . . . , g_p ∈ \mathbb{R}^n$
such that

```math
Z = \{ c + \sum_{i=1}^{p} \xi_i g_i | \xi_i ∈ [−1, 1]\}.
```
It is common to note $Z = (c, \langle g_1 . . . , g_p \rangle)$ or simply
$Z = (c, G)$, where $g_i$ is the $i$-th column of $G$.


 We introduce the order of a zonotope as $o = \frac{k}{n}$.
The order of a Zonotope `Z` in LazySets can be calculated using the function `order(Z)`.

```@example zonotope_example_1
using LazySets, Plots
Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])
plot(Z)
quiver!(fill(1., 4), fill(1., 4), quiver=(genmat(Z)[1, :], genmat(Z)[2, :]), color=:black)
```

There are other useful characterization of zonotopes. A zonotope can be seen as the Minkowski addition of line segments resulting in centrally symmetric convex polytopes as shown in the following figure, which illustrates how each generator spans the zonotope.

The cost is measured in terms of the number of binary operations, $\mathrm{Op}(\cdot)$.


$Z_1 = (c, \langle v_1, \dotsb, v_k \rangle), Z_2 = (d, \langle w_1, \dotsb, w_m \rangle) \subset \mathbb{R}^n, M \in \mathbb{R}^{m \times n}$

$Z_1 \oplus Z_2 = (c+d, \langle v_1, \dotsb, v_k, w_1, \dotsb, w_m \rangle)$

$MZ_1 = (Mc, \langle Mv_1, \dotsb, Mv_k \rangle)$

$CH(Z_1, e^{A\delta}Z_1) \subseteq \frac{1}{2}(c + e^{A\delta}c,\langle v_1 + e^{A\delta}v_1, \dotsb, v_k+e^{A\delta}v_k, v_1 - e^{A\delta}v_1, v_k - e^{A\delta}v_k, c - e^{A\delta}c \rangle )$

| Operation                 | Simplification Rule | Cost               |
|---------------------------|---------------------|--------------------|
| $Z_1 \oplus Z_2$          |                     |      $n$            |
| $MZ_1$                    |                     |      $2mn(k+1)$              |
| $CH(Z_1, e^{A\delta}Z_1)$ |                     |      $2n^2(k+1)+2n(k+2)$         |



### Hausdorff distance

```math
  d_H(\mathcal{X}, \mathcal{Y}) = \max \left( \sup_{x \in \mathcal{X}}\inf_{y \in \mathcal{Y}} \Vert x - y \Vert, \sup_{y \in \mathcal{Y}}\inf_{x \in \mathcal{X}} \Vert x - y \Vert \right)
```


### Taylor models



```@example
#=
using ReachabilityAnalysis, Plots

f(x) = -6x^3 + (13/3)x^2 + (31/3)x
dom = -3.5 .. 3.5

plot(f, -3.5, 3.5, lab="f", xlab="x")

x = Taylor1(5)
set_taylor1_varname("x")
f(x)

rem = 0 .. 0
x0 = 0.0
dom = -3.5 .. 3.5
tm = TaylorModel1(f(x), rem, x0, dom)
=#
```

## Reach-sets

We define the reachable set associated to a time interval $[0, δ]$,
also known as the *flowpipe*, as
```math
\mathcal{F}([0, δ]) = ⋃_{t \in [0, δ]} \mathcal{R}(t).
```
Reachability methods are used to compute rigorous approximations of the flowpipe
for continuous or hybrid systems, in bounded time or unbounded time horizon.
Here we use the term *rigorous* in the formal, or mathematical sense, that no
solution "escapes" the flowpipe, for any trajectory that satisfies the constraints
(initial states, inputs, and noise).

On the other hand, the amount of computation required depends heavily on the
particular problem statement. One notable example is *safety verification*,
which simply stated requires to prove that the flowpipe does not intersect a region
of "bad states". In this setting, one can often reason about the flowpipe lazily,
i.e. without actually computing it in full.


We consider as a running example in this section the simple harmonic oscillator,

## Flowpipes

A flowpipe represents a collection of reach-sets and behaves like their set union.

## Approximations

Given a set ``\mathcal{X} \subseteq \mathbb{R}^n``, any subset of ``\mathcal{X}``
is said to be an *underapproximation*. Conversely, any set containing ``\mathcal{X}``
is said to be an *overapproximation*. In this section we recall the definitions
and give some examples of the basic set operations commonly used to construct
reachability algorithms. Such operations are not only required to propagate
reachable sets

We begin with the notion of Hausdorff distance which, as a way to *measure*
(in the informal sense) the distance between sets, constitutes a practical tool to
quantify the quality of an approximation.


## References
