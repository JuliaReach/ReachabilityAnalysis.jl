```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

## Hausdorff distance

The notion of Hausdorff distance can be used to *measure* the distance between sets.
It constitutes a practical theoretical tool to quantify the quality of an approximation.

```math
  d_H(\mathcal{X}, \mathcal{Y}) = \max \left( \sup_{x ∈ \mathcal{X}}\inf_{y ∈ \mathcal{Y}} \Vert x - y \Vert, \sup_{y ∈ \mathcal{Y}}\inf_{x ∈ \mathcal{X}} \Vert x - y \Vert \right)
```
