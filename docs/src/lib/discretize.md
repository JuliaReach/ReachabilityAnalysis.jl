```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# Discretization

```@contents
Pages = ["discretize.md"]
Depth = 3
```

## Discretize API

```@docs
discretize
```

## Approximation models

```@docs
ReachabilityAnalysis.AbstractApproximationModel
Forward
Backward
CorrectionHull
NoBloating
```

## Exponentiation

```@docs
ReachabilityAnalysis._exp
ReachabilityAnalysis.Φ₁
ReachabilityAnalysis.Φ₂
```
