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

B = BallInf(zeros(2), 1.0) # unit cube
M = [0 1; -1 0] # an affine transformation

Z = linear_map(M, B) # zonotope
```

Zonotopes are commonly represented using their generator representation.
Here, a zonotope $Z ⊆ \mathbb{R}^n$ is defined by a center $c ∈ \mathbb{R}^n$ and
a finite number of generators $g_1, . . . , g_p ∈ \mathbb{R}^n$ such that

```math
Z = \left\{ c + \sum_{i=1}^{p} \xi_i g_i | \xi_i ∈ [−1, 1]\right\}.
```
It is common to note $Z = (c, \langle g_1 . . . , g_p \rangle)$ or simply
$Z = (c, G)$, where $g_i$ is the $i$-th column of $G$. In the examle of above,

```@example zonotope_definition
@show center(Z)

@show genmat(Z)
```

The generators can be

```@example
Z3 = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
```

The order of a zonotope is the ratio between the number of dimensions and the number
of generators; i.e. $o = \frac{k}{n}$. Use the function `order(Z)` to get the order
of a given zonotope.



```@example zonotope_example_1

Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])
plot(Z)
quiver!(fill(1., 4), fill(1., 4), quiver=(genmat(Z)[1, :], genmat(Z)[2, :]), color=:black)
```

## Other characterizations

There are other useful characterization of zonotopes. A zonotope can be seen as the Minkowski addition of line segments resulting in centrally symmetric convex polytopes as shown in the following figure, which illustrates how each generator spans the zonotope.

## Operations with zonotopes

The cost is measured in terms of the number of binary operations, $\mathrm{Op}(\cdot)$.


$Z_1 = (c, \langle v_1, \dotsb, v_k \rangle), Z_2 = (d, \langle w_1, \dotsb, w_m \rangle) \subset \mathbb{R}^n, M \in \mathbb{R}^{m \times n}$

$Z_1 \oplus Z_2 = (c+d, \langle v_1, \dotsb, v_k, w_1, \dotsb, w_m \rangle)$

$MZ_1 = (Mc, \langle Mv_1, \dotsb, Mv_k \rangle)$

$CH(Z_1, e^{A\delta}Z_1) \subseteq \frac{1}{2}(c + e^{A\delta}c,\langle v_1 + e^{A\delta}v_1, \dotsb, v_k+e^{A\delta}v_k, v_1 - e^{A\delta}v_1, v_k - e^{A\delta}v_k, c - e^{A\delta}c \rangle )$

| Operation                 | Simplification Rule | Cost               |
|---------------------------|---------------------|--------------------|
| $Z_1 \oplus Z_2$          |                     |      $n$            |
| $MZ_1$                    |                     |      $2mn(k+1)$              |
| $CH(Z_1, e^{A\delta}Z_1)$ |                     |      $2n^2(k+1)+2n(k+2)$         |
