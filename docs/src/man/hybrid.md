```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Hybrid systems

## Introduction


Our running example is the *bouncing ball* model; although it is a very hybrid automaton,
it can be used to introduce the main notions involved in hybrid systems reachability.

![Hybrid automaton of the bouncing ball model](../assets/bouncing_ball_annotations.png)

## Formalism

TODO: definition of HA

## Hybrid solver algorithm

TODO: high-level description of the hybrid solve

## Flowpipe-guard intersections

In this section we illustrate the available methods to compute flowpipe-guard intersections. We will see how to use the function `cluster` to select a clustering strategies to cope with the case when there are several intersecting sets. We will also show examples of choosing different intersection templates.

We'll take for our running example a two-dimensional rotating system with dynamics

```math
 \begin{aligned}
   \dot{x} &= y \\
   \dot{y} &= - x
\end{aligned}
```
and define the guard ``G: \\{(x, y) \in \\mathbb{R}^2: x ≥ 1.3 \\}``.

## Customizing the solver

TODO: further examples on how to customize the discrete post-operator
