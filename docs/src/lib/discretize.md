```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Discretization

```@contents
Pages = ["discretize.md"]
Depth = 3
```

## Discretize API

```@docs
normalize
discretize
```

## Approximation models

```@docs
NoBloating
Forward
Backward
FirstOrder
FirstOrderZonotope
SecondOrderddt
CorrectionHull
StepIntersect
```

## Exponentiation

The state transition matrix of the linear ODE $x'(t) = Ax(t) + u(t)$ at time
$\delta > 0$ is $\Phi = e^{A\delta}$, hence the algorithms usually require to
compute exponential matrices. There are distinct ways to compute the matrix
exponential $e^{A\delta}$ depending on the type of $A$
(see e.g. [^HIH08]). The available methods can be used through the (unexported) function `_exp`.

For high dimensional systems (typically `n > 2000`), computing the matrix exponential
is expensive hence it is preferable to compute the action of the matrix exponential
over vectors when needed, that is, $e^{δA} v$ for each $v$. This method is particularly
well-suited if `A` is vert sparse. Use the option `exp=:krylov` (or `exp=:lazy`) for this purpose.

```@docs
ReachabilityAnalysis.Exponentiation
ReachabilityAnalysis.Exponentiation._exp
ReachabilityAnalysis.Exponentiation.Φ₁
ReachabilityAnalysis.Exponentiation.Φ₂
ReachabilityAnalysis.Exponentiation.AbstractExpAlg
ReachabilityAnalysis.Exponentiation.BaseExpAlg
ReachabilityAnalysis.Exponentiation.LazyExpAlg
ReachabilityAnalysis.Exponentiation.IntervalExpAlg
ReachabilityAnalysis.Exponentiation.PadeExpAlg
```

## References

[^HIH08]: Higham, Nicholas J. Functions of matrices: theory and computation. Society for Industrial and Applied Mathematics, 2008.
