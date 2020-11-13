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

### Representing sets with support functions

Support functions are one of the central tools in set-based reachability, because

 Let $\mathcal{X}$ be a non-empty subset of a finite dimensional space $\mathcal{X} \subseteq \mathbb{R}^n$. The *support function* of $\mathcal{X} \subseteq \mathbb{R}^n$ attributes to each vector $d \in \mathbb{R}^n$ the real number

```math
  ρ(d, \mathcal{X}) = \max \{ d^T x : x \in \mathcal{X} \},
```
where $d^T$ denotes the transpose of the (column) vector $d$.

TODO: example with triangle
A convex set can be represe

The

### Polyhedra

### Zonotopes

Zonotopes are a sub-class of polytopes defined as the image of a unit cube under


### Hausdorff distance



### LazySets


### Taylor models


## Reach-sets



## Flowpipes

A flowpipe represents a collection of reach-sets that behaves like a set union.
