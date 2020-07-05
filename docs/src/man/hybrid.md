# Hybrid systems

## Introduction


Our running example is the *bouncing ball* model; although it is a very hybrid automaton, it can be used to introduce the main notions involved in hybrid systems reachability.


## Clocked linear dynamics

So far we have focused on transitions that involve "spatial" variables.
If the system under consideration has transitions governed by time variables,
i.e. by variables whose dynamics are of the form ``t' = 1``, then decoupling the
spatial variables with the clock variables gives a computational advantage.
We refer to [[HG19]].

## Flowpipe-guard intersection

In this section we illustrate the available methods to compute flowpipe-guard intersections. We will see how to use the function `cluster` to select a clustering strategies to cope with the case when there are several intersecting sets. We will also show examples of choosing different intersection templates.

We'll take for our running example a two-dimensional rotating system with dynamics

```math
 \begin{aligned}
   \dot{x} &= y \\
   \dot{y} &= - x
\end{aligned}
```
and define the guard ``G: \\{(x, y) \in \\mathbb{R}^2: x ≥ 1.3 \\}``.
