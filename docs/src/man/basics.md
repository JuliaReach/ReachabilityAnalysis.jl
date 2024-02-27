```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

# Set representations

## Introduction

In simple terms, reachability analysis is concerned with studying the sets of states
that a system can reach, starting from a set of initial states and under the influence
of a set of input trajectories and parameter values.

To fix ideas, consider a purely continuous dynamical system over some state
space ``\mathcal{X} ⊆ \mathbb{R}^n`` defined via a differential equation of the form:

```math
    x'(t) = f(x, u),
```
where ``u(t) ∈ \mathcal{U}(t)`` ranges over some specified set of admissible input signals.
Given a set of initial states ``\mathcal{X}_0 ⊆ \mathcal{X}``, the problem we are
interested in is to compute *all* the states visited by trajectories of the system starting
from any ``x_0 ∈ \mathcal{X}_0`` and for any admissible input function ``u(t)``. In other
sections of this manual we also consider generalizations of this problem such as uncertain
parameters, or dynamical systems which are hybrid, i.e. mixing continuous dynamics and discrete transitions.

In order to compute with sets, we need to find convenient ways to represent and operate with them.
In the next section we consider some common set representations and show how to work with `LazySets.jl`.


----

More technically, we define the reachable set at a given time point
``δ ∈ \mathbb{R}``, also known as the *reach-set* for the ODE
```math
x'(t) = f(x(t), u(t), p(t), t),
```
as given by
```math
\mathcal{R}(δ) := \left\{ x(δ) = ∫_0^δ f(x(t), u(t), p(t), t) dt, x(0) ∈ X_0, u(t) ∈ \mathcal{U}, p(t) ∈ \mathcal{P} \right\}.
```
Here ``X_0`` denotes the set of initial states, ``\mathcal{U}`` denotes the input set,
and ``\mathcal{P}`` denotes the parameter values. For practical problems, the set
``\mathcal{R}(δ)`` cannot be obtained exactly, and reachability methods aim at
computing suitable over-approximations (or under-approximations) of it.


## Set representations

Subsets of ``\mathbb{R}^n`` can be represented in different ways. Depending on the
type of operation that we want to apply, one set representation may be more convenient
than another one. There are two important characteristics:

- **Closure:** A set representation is said to be *closed* under a given operation if the result of the operation is again a set in the same class.

- **Cost:** To measure the cost of making a given operation on a set, we consider the total number of binary operations, denoted with ``\mathrm{Op}(⋅)``. (....)

In the rest of this section we define the usual set representations used in reachability analysis and recall some fundamental properties. Moreover, we show how to define these sets using the library `LazySets.jl`.

## Over- and under- approximations


Given a set ``\mathcal{X} ⊆ \mathbb{R}^n``, any subset of ``\mathcal{X}``
is said to be an *underapproximation*. Conversely, any set containing ``\mathcal{X}``
is said to be an *overapproximation*.

### Polyhedra

A hyperplane is the set ``\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx = b\}``, where ``a ∈ \mathbb{R}^n`` is the normal vector and ``b ∈ \mathbb{R}`` is the displacement.

In `LazySets.jl`, the type `hyperplane` can be used to define a hyperplane.


A halfspace is the set ``\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx ≤ b\}``, where ``a ∈ \mathbb{R}^n`` is the normal vector and ``b ∈ \mathbb{R}`` is the displacement. Please note that a half-space defines the region on one side of the hyperplane ``a^Tx = b``.

In `LazySets.jl`, the type `HalfSpace` can be used to define a half-space.


When we consider the intersection of a finite subset of half-spaces, we get a polyhedron. A polyhedron is thus the set defined as ``\mathcal{P} ⊆ \mathbb{R}^n``,
```math
    \mathcal{P} = \{x ∈ \mathbb{R}^n | \bigcap_{i=1}^m a^T_i x ≤ b_i \},
```
where ``a_i ∈ \mathbb{R}^n`` and ``b_i ∈ \mathbb{R}``. A polytope is a bounded polyhedron.

In `LazySets.jl`, the types `HPolyhedron` and `HPolytope` representent polyhedron and polytopes respectively. For example, . . . .


### Support functions

The support function of a compact set ``\mathcal{X} ⊆ \mathbb{R}^n`` attributes to a direction ``ℓ ∈ \mathbb{R}^n`` the real number
```math
    ρ_{\mathcal{X}} (ℓ) = \max\{ ℓ^T x | x ∈ \mathcal{X} \}.
```
It is important to note that for a given direction ``ℓ``, the support function defines the position of a halfspace
```math
    \mathcal{H}_{ℓ} = \{x ∈ \mathbb{R}^n | ℓ^T x ≤ ρ_{\mathcal{X}}(ℓ)\}.
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

- ``ρ_{\mathcal{X} \oplus \mathcal{Y}}(ℓ) =

- ``ρ_{M\mathcal{X}}(ℓ) =  for any
Here ``M\mathcal{X}`` represents the linear map for ``M ∈ \mathbb{R}^{m \times n}``.
- ``ρ_{CH(\mathcal{X}, \mathcal{Y})}(ℓ) =

The following table summarizes the number of
| Operation                                     | Simplification rule | Cost  |
|-----------------------------------------------|-------|------|
| ``ρ_{\mathcal{X} \oplus \mathcal{Y}}(ℓ)`` |``ρ_{\mathcal{X}}(ℓ) + ρ_{\mathcal{Y}}(ℓ)``  |  ``1``     |
| ``ρ_{M\mathcal{X}}(ℓ)``                   |``ρ_{\mathcal{X}}(M^Tℓ)``|   ``2mn``   |
| ``ρ_{CH(\mathcal{X}, \mathcal{Y})}(ℓ)``   |``\max\{ρ_{\mathcal{X}(ℓ)}, ρ_{\mathcal{Y}(ℓ)}\}``|  ``1``    |

which touches and contains ``\mathcal{X}``. If ``ℓ`` is of unit length, then
``ρ_{\mathcal{X}}(ℓ)`` is the signed distance of ``\mathcal{H}_{ℓ}`` to the origin.

Evaluating the support function for a set of directions ``L ⊆ \mathbb{R}^n`` provides an overapproximation

```math
    ⌈ \mathcal{X} ⌉ _L = \bigcap_{ℓ ∈ L} \{ x ∈ \mathbb{R}^n | ℓ^T x ≤ ρ_{\mathcal{X}}(ℓ) \}
```

i.e., ``\mathcal{X} ⊆ ⌈ \mathcal{X} ⌉ _L``. If ``L = \mathbb{R}^n``, then ``\mathcal{X} = ⌈ \mathcal{X} ⌉ _L``, so the support function represents any convex set ``\mathcal{X}`` exactly. If ``L`` is a finite set of directions ``L = {ℓ_1, . . . , ℓ_m}``, then ``⌈ \mathcal{X} ⌉ _L`` is a polyhedron.






## Reach-sets

We define the reachable set associated to a time interval ``[0, δ]``,
also known as the *flowpipe*, as
```math
\mathcal{F}([0, δ]) = ⋃_{t ∈ [0, δ]} \mathcal{R}(t).
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
