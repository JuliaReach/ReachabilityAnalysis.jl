# Systems

Systems types are defined in the library [MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl). Apart from purely discrete or continuous, hybrid automata
(i.e. those mixing discrete-continuous dynamics) are defined in
[HybridSystems.jl](https://github.com/blegat/HybridSystems.jl).

```@contents
Pages = ["systems.md"]
Depth = 3
```

## Types and macros

The API reference for systems types and macros can be found in the
[MathematicalSystems.jl](https://juliareach.github.io/MathematicalSystems.jl/latest/man/systems/)
documentation. Two commonly used macros are `@system` and `@ivp`, used to
define a system and an initial-value problem respectively.

## Normalization

```@docs
normalize(::AbstractSystem)
normalize(::InitialValueProblem)
ReachabilityAnalysis.add_dimension
```

## Homogeneization

```@docs
homogenize
```

## Hybrid systems

```@docs
HACLD1
DiscreteTransition
constrained_dimensions(::HybridSystem)
```
