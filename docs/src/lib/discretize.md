```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# Discretization

```@contents
Pages = ["discretize.md"]
Depth = 3
```

## Approximation models

```@docs
ReachabilityAnalysis.AbstractApproximationModel
Forward
Backward
CorrectionHull
NoBloating
```

## Discretize API

```@docs
discretize
```

## Exponentiation

```@docs
ReachabilityAnalysis._exp
ReachabilityAnalysis.Φ₁
ReachabilityAnalysis.Φ₂
```
