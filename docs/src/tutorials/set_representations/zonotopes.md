```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

## Zonotopes

Zonotopes are a sub-class of polytopes defined as the image of a unit cube under
an affine transformation. For example,

```@example zonotope_definition
using LazySets, Plots #hide
using LazySets: center #hide

B = BallInf(zeros(2), 1.0) # unit cube
M = [0 1; -1 0] # an affine transformation

Z = linear_map(M, B) # zonotope
```

Zonotopes are commonly represented using their generator representation.
Here, a zonotope ``Z ⊆ \mathbb{R}^n`` is defined by a center ``c ∈ \mathbb{R}^n`` and
a finite number of generators ``g_1, . . . , g_p ∈ \mathbb{R}^n`` such that

```math
Z = \left\{ c + \sum_{i=1}^{p} ξ_i g_i | ξ_i ∈ [−1, 1]\right\}.
```
It is common to note ``Z = (c, ⟨ g_1 . . . , g_p ⟩)`` or simply
``Z = (c, G)``, where ``g_i`` is the ``i``-th column of ``G``. In the example above,

```@example zonotope_definition
@show center(Z)

@show genmat(Z)
```

The generators can be

```@example zonotope_definition
Z3 = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
```

The order of a zonotope is the ratio between the number of dimensions and the number
of generators; i.e. ``o = \frac{k}{n}``. Use the function `order(Z)` to get the order
of a given zonotope.



```@example zonotope_definition
Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])
plot(Z)
quiver!(fill(1., 4), fill(1., 4), quiver=(genmat(Z)[1, :], genmat(Z)[2, :]), color=:black)
```

## Other characterizations

There are other useful characterization of zonotopes. A zonotope can be seen as the Minkowski addition of line segments resulting in centrally symmetric convex polytopes as shown in the following figure, which illustrates how each generator spans the zonotope.

## Operations with zonotopes

The cost is measured in terms of the number of binary operations, ``\mathrm{Op}(⋅)``.


``Z_1 = (c, ⟨ v_1, …, v_k ⟩), Z_2 = (d, ⟨ w_1, …, w_m ⟩) ⊆ \mathbb{R}^n, M ∈ \mathbb{R}^{m \times n}``

``Z_1 \oplus Z_2 = (c+d, ⟨ v_1, …, v_k, w_1, …, w_m ⟩)``

``MZ_1 = (Mc, ⟨ Mv_1, …, Mv_k ⟩)``

``CH(Z_1, e^{Aδ}Z_1) ⊆ \frac{1}{2}(c + e^{Aδ}c,⟨ v_1 + e^{Aδ}v_1, …, v_k+e^{Aδ}v_k, v_1 - e^{Aδ}v_1, v_k - e^{Aδ}v_k, c - e^{Aδ}c ⟩ )``

| Operation              | Simplification Rule | Cost                  |
|------------------------|---------------------|-----------------------|
| ``Z_1 \oplus Z_2``     |                     | ``n``                 |
| ``MZ_1``               |                     | ``2mn(k+1)``          |
| ``CH(Z_1, e^{Aδ}Z_1)`` |                     | ``2n^2(k+1)+2n(k+2)`` |
