```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

# Basics

## Problem statement

In simple terms, reachability analysis is concerned with studying the sets of states
that a system can reach, starting from a set of initial states and under the
influence of a set of input trajectories and parameter values.
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

Up to now we have discussed about the continuous case only, but there is a rich
literature in hybrid systems reachability; *hybrid* here means those dynamical
systems which are given by one or more continuous-time dynamics (often, systems
of ODEs in each mode or location) coupled with discrete transitions between
continuous modes. In our context it is standard to model these systems using the
terminology of *hybrid automata*, and we also model hybrid systems with such framework
in this library. The concept of reach-set, flowpipe and safety verification are
naturally extended to hybrid automata, although there is the additional complication
that the flowpipe must include the behaviors for all possible transitions between
discrete modes that are compatible with the dynamics.

## Set operations

Given a set ``\mathcal{X} \subseteq \mathbb{R}^n``, any subset of ``\mathcal{X}``
is said to be an *underapproximation*. Conversely, any set containing ``\mathcal{X}``
is said to be an *overapproximation*. In this section we recall the definitions
and give some examples of the basic set operations commonly used to construct
reachability algorithms. Such operations are not only required to propagate
reachable sets

We begin with the
notion of Hausdorff distance which, as a way to *measure* (in the informal sense)
the distance between sets, constitutes a practical tool to quantify the quality
of an approximation.

### Hausdorff distance

```math
  d_H(\mathcal{X}, \mathcal{Y}) = \max \left( \sup_{x \in \mathcal{X}}\inf_{y \in \mathcal{Y}} \Vert x - y \Vert, \sup_{y \in \mathcal{Y}}\inf_{x \in \mathcal{X}} \Vert x - y \Vert \right)
```

## Set representations

In this section we consider

### Zonotope set representation

A zonotope $Z ⊆ \mathbb{R}^n$ is defined by a center $c ∈ \mathbb{R}^n$ and a finite number of generators $v1, . . . , vk ∈ \mathbb{R}^n$:
```math
Z = \{ c + \sum_{k}^{i=1} α_i v_i | α_i ∈ [−1, 1]\}.
```
A common denotation for zonotopes is $Z = (c, \langle v1, . . . , vk \rangle)$. We introduce the order of a zonotope as $o = \frac{k}{n}$.
The order of a Zonotope `Z` in LazySets can be calculated using the function `order(Z)`.
A zonotope can be seen as the Minkowski addition of line segments resulting in centrally symmetric
convex polytopes as shown in the following figure, which illustrates how each generator spans the zonotope.

$Z_1 = (c, \langle v_1, \dotsb, v_k \rangle), Z_2 = (d, \langle w_1, \dotsb, w_m \rangle) \subset \mathbb{R}^n, M \in \mathbb{R}^{m \times n}$

$Z_1 \oplus Z_2 = (c+d, \langle v_1, \dotsb, v_k, w_1, \dotsb, w_m \rangle)$

$MZ_1 = (Mc, \langle Mv_1, \dotsb, Mv_k \rangle)$

$CH(Z_1, e^{A\delta}Z_1) \subseteq \frac{1}{2}(c + e^{A\delta}c,\langle v_1 + e^{A\delta}v_1, \dotsb, v_k+e^{A\delta}v_k, v_1 - e^{A\delta}v_1, v_k - e^{A\delta}v_k, c - e^{A\delta}c \rangle )$

| Operation                 | Cost                |
|---------------------------|---------------------|
| $Z_1 \oplus Z_2$          | $n$                 |
| $MZ_1$                    | $2mn(k+1)$          |
| $CH(Z_1, e^{A\delta}Z_1)$ | $2n^2(k+1)+2n(k+2)$ |

```@example zonotope_example_1
using LazySets, Plots
Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])
plot(Z)
quiver!(fill(1., 4), fill(1., 4), quiver=(genmat(Z)[1, :], genmat(Z)[2, :]), color=:black)
```

### Representing sets with support functions


In this subsection, we provide definitions for polyhedra and support functions, and recall some fundamental properties.
A halfspace $\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx ≤ b\}$, with normal vector $a ∈ \mathbb{R}^n$ and $b ∈ \mathbb{R}$ is one half of the space after dividing it by a hyperplane.
A polyhedron $\mathcal{P} ⊆ \mathbb{R}^n$ is the intersection of a finite number of halfspaces, written as
```math
    \mathcal{P} = \{x \in \mathbb{R}^n | \bigcap_{i=1}^m a^T_i x \leq b_i \},
```
where $a_i \in \mathbb{R}^n$ and $b_i \in \mathbb{R}$. A polytope is a bounded polyhedron. The support function of a compact set $\mathcal{X}$ attributes to a direction $\ell \in \mathbb{R}^n$
```math
    \rho_{\mathcal{X}} (\ell) = max\{ \ell^T x | x \in \mathcal{X} \}.
```
for a given direction $\ell$, it defines the position of a halfspace
```math
    \mathcal{H}_{\ell} = \{x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell)\},
```




$\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell) = \rho_{\mathcal{X}}(\ell) + \rho_{\mathcal{Y}}(\ell)$

$\rho_{M\mathcal{X}}(\ell) = \rho_{\mathcal{X}}(M^T\ell), (M \in \mathbb{R}^{m \times n})$

$\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell) = max{\rho_{\mathcal{X}(\ell)}, \rho_{\mathcal{Y}(\ell)}}$

| Operation                                     | Cost  |
|-----------------------------------------------|-------|
| $\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell)$ | $1$   |
| $\rho_{M\mathcal{X}}(\ell)$                   | $2mn$ |
| $\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell)$   | $1$   |

which touches and contains $\mathcal{X}$ . If $\ell$ is of unit length, then
$\rho_{\mathcal{X}}(\ell)$ is the signed distance of $\mathcal{H}_{\ell}$ to the origin.
Evaluating the support function for a set of directions $L ⊆ \mathbb{R}^n$ provides an overapproximation
```math
    \lceil \mathcal{X} \rceil _L = \bigcap_{\ell \in L} \{ x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell) \}
```

i.e., $\mathcal{X} ⊆ \lceil \mathcal{X} \rceil _L$. If $L = \mathbb{R}^n$, then $\mathcal{X} = \lceil \mathcal{X} \rceil _L$, so the support function represents any convex set $\mathcal{X}$ exactly. If $L$ is a finite set of directions $L = {\ell_1, . . . , \ell_m}$, then $\lceil \mathcal{X} \rceil _L$ is a polyhedron.

### Polyhedra

### Zonotopes

Zonotopes are a sub-class of polytopes defined as the image of a unit cube under


### Hausdorff distance



### LazySets


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

We consider as a running example in this section the simple harmonic oscillator,

## Flowpipes

A flowpipe represents a collection of reach-sets that behaves like a set union.
