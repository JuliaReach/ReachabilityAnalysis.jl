## BFFPSV18

### Reachable states approximation

In a nutshell, we overapproximate the reachable states of an affine system by
solving a set-based recurrence. The key idea is that we first decompose the system into (low-dimensional)
subsystems and later compose the results as a Cartesian product.
Thus we have to solve many cheap problems instead of one hard problem.
Since solving the recurrence scales superlinearly with the dimension, this
approach is very scalable.

### Decomposition error

Decomposition typically involves a loss in precision, and so does this approach.
The good thing is that we can decompose the recurrence as well, which allows us
to analyze each of the subsystems independently by only referring to the initial
states of the other subsystems.
Consequently, there are two main sources for precision loss:
1. Decomposition of the initial states: If two subsystems are interdependent initially.
2. Representation of the reachable states as a Cartesian product: If two subsystems are interdependent in the dynamics.
3. Representation of the reachable states in general: The reachable states of affine systems cannot be represented precisely in all cases. This is a problem that all approaches suffer from. We overapproximate the reachable states by (unions of) convex polytopes.

### Checking safety properties

The problem of checking a safety property can be reduced to a reachability
problem.
We provide special support for this reduction by inlining the property check
into the reachable states computation.
This has two benefits:
1. We fail fast when the property is violated in our abstraction.
2. The check is usually cheaper than computing the full reachable states. This is because we are often only interested in an upper or lower bound of a variable.

## Lazy sets

To represent sets of states, we use the `LazySets` package which provides
exact but lazy (i.e. symbolic) representations of common sets.
