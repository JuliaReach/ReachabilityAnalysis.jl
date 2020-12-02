```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Exploiting structure

In this section we consider reachability analysis for subsets of variables
of a given linear system. Let

```math
x'(t) = Ax(t) + v(t)
```

```@example
a = 1
b = 2
a + b
```

IDEAS:

- Method to lazily compute the flowpipe when we are only interested in outputs. Example with linear combination of state varibles with LGG09 (eg. some of SLICOT benchmaris properties).
- Show case using Krylov when the state matrix is sparse.
- Show option when Phi is sparse.
- A section about BFFPSV18.
