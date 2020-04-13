```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# Discretization

```@contents
Pages = ["discretize.md"]
Depth = 3
```

## Linear systems

TODO: document ReachabilityAnalysis.discretize

## Exponentiation


```@docs
ReachabilityAnalysis._exp
ReachabilityAnalysis.Φ₁
ReachabilityAnalysis.Φ₂
```

Let $A ∈ \mathbb{R}^{n×n}$ and for $t ≥ 0$ consider the integral
$\int_0^t e^{Aξ}dξ$. If $A$ is invertible, this integral (...)
