```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# Discretization

```@contents
Pages = ["discretize.md"]
Depth = 3
```

## Linear systems

```@docs
ReachabilityAnalysis.discretize
ReachabilityAnalysis.discretize_interpolation
ReachabilityAnalysis.discretize_firstorder
ReachabilityAnalysis.discretize_nobloating
ReachabilityAnalysis.discretize_interval_matrix
```

## Exponentiation

Let $A ∈ \mathbb{R}^{n×n}$ and for $t ≥ 0$ consider the integral
$\int_0^t e^{Aξ}dξ$. If $A$ is invertible, this integral


```@docs
ReachabilityAnalysis.exp_Aδ
ReachabilityAnalysis.Φ₁
ReachabilityAnalysis.Φ₂
```
