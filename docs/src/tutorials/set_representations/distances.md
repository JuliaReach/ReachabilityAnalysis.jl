```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

## Hausdorff distance

Tthe notion of Hausdorff distance can be used to *measure* the distance between sets.
It constitutes a practical theoretical tool to quantify the quality of an approximation.

```math
  d_H(\mathcal{X}, \mathcal{Y}) = \max \left( \sup_{x \in \mathcal{X}}\inf_{y \in \mathcal{Y}} \Vert x - y \Vert, \sup_{y \in \mathcal{Y}}\inf_{x \in \mathcal{X}} \Vert x - y \Vert \right)
```
