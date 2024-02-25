```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

## Modeling uncertain inputs

In this section of the manual we consider systems with linear continuous dynamics

```math
x'(t) = Ax(t) + Bu(t),\qquad x(0) ∈ X_0,\qquad (1)
```
where ``A ∈ \mathbb{R}^{n\times n}``, ``x(t) ∈ X ⊆ \mathbb{R}^n``,
``B ∈ \mathbb{R}^{n\times m}`` and ``u(t) ∈ U ∈ \mathbb{R}^m``.
Systems with mixed discrete/continuous dynamics of type (1) are discussed
in the section [hybrid_systems](@ref).

Generally, we distinguish between three types of inputs:

1. Fixed inputs, where ``u(t)`` is precisely known. In some cases, ``u(t)`` is a constant;
   in other cases, it is a function specified by the user.
