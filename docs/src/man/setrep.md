```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
    using Plots
end
CurrentModule = ReachabilityAnalysis
```

# Set representations

Computational problems that involve subsets of $\mathbb{R}^n$ can be formulated
and solved in many different ways, and each approach has its advantages and
disadvantages: there is no "one size fits all" set represententation in our discipline.
Depending on the question, one set representation may be more convenient than another one.
In this section we define the most common set representations used in reachability
analysis and recall some fundamental properties. Various examples are presented using
[LazySets.jl](https://github.com/JuliaReach/LazySets.jl), our core package for set representations.

Algorithms that use [LazySets.jl](https://github.com/JuliaReach/LazySets.jl) can
explore different approaches with minimal changes in their code, because conversion
between set representations as well as overapproximation or underapproximation algorithms
are implemented. Such transformations are discussed in this end of this section.

## Lazy and concrete representations

[LazySets.jl](https://github.com/JuliaReach/LazySets.jl) is our core library for
set-based computations in the Euclidean space. A few dozen set representations are defined, such as
polyhedra in constraint and in vertex representation, ellipsoids, balls in different norms,
specific classes of polyhedra such as intervals, zonotopes or hyperrectangles, for which
specialized methods apply. The library can handle the usual operations between sets,
such as Minkowski sum, cartesian product, convex hull, intersection and unions, to mention
only a few. There are two key aspects of the library: substantial algorithm specialization by
using Julia's multiple-dispatch paradigm; and an extensive use of the
[lazy evaluation](https://en.wikipedia.org/wiki/Lazy_evaluation) paradigm:
set computations can be done either lazily or eagerly (or a combiation of them),
but by design, all operations are lazy by default, i.e.
`X ⊕ MY` returns the *lazy* Minkowski sum of `X` and the *lazy* linear map of `Y`
under `M`. Hence, formulating the computation `Z = X ⊕ MY` has "cost zero" until
we want to obtain e.g. the support function of `Z` along a given direction, or we
want to convert it to a concrete representation, e.g. a zonotope. We discuss these
points in the next paragraphs.

As mentioned, of the key features of [LazySets.jl](https://github.com/JuliaReach/LazySets.jl) is specialization.
On a more technical level, the library is highly optimized to get the best possible performance given the following
two restrictions: the set type and the set dimension. For instance, consider the linear map
transformation, that is, to compute $Y = \{y : y = Ax for some x ∈ X\}$.
In [LazySets.jl](https://github.com/JuliaReach/LazySets.jl), `linear_map(A, X)` returns a set representation of Y and the algorithm that is actually used, as well as the type of `Y`, depends on the types
of its arguments and the dimension of `X`, that is, multiple dispatch.
For example, if `X` is a polygon in vertex representation (`VPolygon`) then the
linear map `A` is applied to each vertex. However, if `X` is a 30-dimensional
polyhedron in halfspace representation (`HPolyhedron`), then the halfspace representation
is used (if possible). In [LazySets.jl](https://github.com/JuliaReach/LazySets.jl),
set operations can be performed in two possible (complementary) ways: concretely (or eagerly)
or lazily. Concrete computation means that the result is a set, and in the lazy computation the
result is an object that represents the operation between the given sets. For instance,
while `linear_map(A, X)` computes the concrete linear map as explained in the previous paragraph,
`LinearMap(A, B)`, or simply A * X, computes the lazy linear map, which is just a new object
which "holds" the linear map computation until it is actually needed. In other words,
`LinearMap(A, X)` can be used to reason about the linear map even if `X` is high-dimensional,
since this command just builds an object representing the linear map of `X` by `A`.
In Julia, Unicode symbols such as `A * X`, `X ⊕ Y`, `X ⊖ Y`, `X × Y`, all default to
lazy operations by design. Unicode is input using the latex name of the symbol
followed by the TAB key, such as `\oplus[TAB]` for the (lazy) Minkowski sum `⊕`.

Finally, one can combine lazy set operations to build lazy expressions that represent
several operations between sets, such as eg. `Q = (Z ⊕ A*X) × T`. By means of the basic
tools of convex geometry, such as support function theory, useful information about
`Q` can be obtained without actually computing the linear map, Minkowski sum and cartesian product
in the above computation. For example, if the computation only involves querying
information about `Q` in a small number of directions, the lazy approach can be done very efficiently.

## Performance measures

When comparing representations for a given task, there are important properties
to be considered, among them the main characteristics are:

- **Convexity:** Some representations can only be used to represent convex sets,
  e.g. zonotopes or template polyhedra. On the other hand, Taylor models or set unions
  can represent non-convex sets.

- **Closure:** A set representation is said to be *closed* under a given operation
  if the result of the operation is again a set in the same class. For example, zonotopes
  are closed under linear map of a zonotope is again a zonotope; on the other hand,
  hyperrectangles are not closed under such operation because the linear map of a hyperrectangle
  is not a hyperrectangle in general (it is a zonotope).

- **Scalability:** To measure the cost of making a given operation on a set in a simplified
  model of comutation, it is common to consider the *time complexity*, associated to
  the total number of binary operations to solve an instance of the computational
  problem as a function of characteristics of the input. A related measure is the
  *space complexity*, associated to the total number of memory space required to solve
  an instance of the computational problem.

- **Dependency preservation:** Set representations are said to be dependency preserving
  under a given set operation if we can reason about subsets of the resulting operation
  by reasoning on the representation domain, i.e. without actually having to
  recompute the operation for the subset.

## Hyperrectangular sets

A special class of polyhedra are (...)

The one-dimensional case are `Interval` representations (...)


```@example hyperrectangle_example_1
H = Hyperrectangle(low=[0.9, -0.1], high=[1.1, 0.1])

plot(H, ratio=1.)
```


```@example hyperrectangle_example_1
typeof(H)
```


## Zonotopic sets

Zonotopes are a sub-class of polytopes defined as the image of a unit cube under
an affine transformation. An equivalent characterization of zonotopes is the
generator representation. Here, $Z ⊆ \mathbb{R}^n$ is defined by a center $c ∈ \mathbb{R}^n$ and a finite number of generators $g_1, . . . , g_p ∈ \mathbb{R}^n$
such that

```math
Z = \{ c + \sum_{i=1}^{p} \xi_i g_i | \xi_i ∈ [−1, 1]\}.
```
It is common to note $Z = (c, \langle g_1 . . . , g_p \rangle)$ or simply
$Z = (c, G)$, where $g_i$ is the $i$-th column of $G$. The order of an $n$-dimensional
zonotope with $p$ generators is $o = \frac{p}{n}$.

For example, we can see a hyperrectangle as special case of a zonotope of order 1,
where the generators are diagonal:

```@example hyperrectangle_example_1
Z = convert(Zonotope, H);

typeof(Z)
```

```@example hyperrectangle_example_1
center(Z)
```


```@example hyperrectangle_example_1
generators(Z)
```



```@example zonotope_example_1
Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])
```

The order of `Z` can be obtained using the function `order(Z)`:

```@example zonotope_example_1
Z = Zonotope([1, 1.], [-1 0.3 1.5 0.3; 0 0.1 -0.3 0.3])

plot(Z)
quiver!(fill(1., 4), fill(1., 4), quiver=(genmat(Z)[1, :], genmat(Z)[2, :]), color=:black)
```

There are other useful characterization of zonotopes. A zonotope can be seen as the Minkowski addition of line segments resulting in centrally symmetric convex polytopes as shown in the following figure, which illustrates how each generator spans the zonotope.

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

## Polyhedra

A hyperplane is the set $\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx = b\}$,
where $a ∈ \mathbb{R}^n$ is the normal vector and $b ∈ \mathbb{R}$ is the displacement.

In `LazySets.jl`, the type `hyperplane` can be used to define a hyperplane. For example, . . . .

polytope is bounded, polyhedron is unbounded (...)

EJEMPLO <<<<<

A halfspace is the set $\mathcal{H} = \{x ∈ \mathbb{R}^n | a^Tx ≤ b\}$, where $a ∈ \mathbb{R}^n$ is the normal vector and $b ∈ \mathbb{R}$ is the displacement. Please note that a half-space defines the region on one side of the hyperplane ``a^Tx = b``.

In `LazySets.jl`, the type `HalfSpace` can be used to define a half-space. For example, . . . .

EJEMPLO <<<<<

When we consider the intersection of a finite subset of half-spaces, we get a polyhedron. A polyhedron is thus the set defined as $\mathcal{P} ⊆ \mathbb{R}^n$,
```math
    \mathcal{P} = \{x \in \mathbb{R}^n | \bigcap_{i=1}^m a^T_i x \leq b_i \},
```
where $a_i \in \mathbb{R}^n$ and $b_i \in \mathbb{R}$. A polytope is a bounded polyhedron.

In `LazySets.jl`, the types `HPolyhedron` and `HPolytope` representent polyhedron and polytopes respectively. For example, . . . .

## Support functions

The support function of a compact set $\mathcal{X} \subseteq \mathbb{R}^n$ attributes to a direction $\ell \in \mathbb{R}^n$ the real number
```math
    \rho_{\mathcal{X}} (\ell) = \max\{ \ell^T x | x \in \mathcal{X} \}.
```
It is important to note that for a given direction $\ell$, the support function defines the position of a halfspace
```math
    \mathcal{H}_{\ell} = \{x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell)\}.
```

```@example
E = Ellipsoid([1.0, 1.0], [2.7 0.28; 0.28 1.04])

H1 = HalfSpace([1.0, 2.0], ρ([1.0, 2.0], E))
H2 = HalfSpace([-1.0, -2.0], ρ([-1.0, -2.0], E))
H3 = HalfSpace([-1.0, 2.0], ρ([-1.0, 2.0], E))
H4 = HalfSpace([1.0, -2.0], ρ([1.0, -2.0], E))

P = HPolytope([H1, H2, H3, H4])

plot(E)
plot!(P)

plot!(Singleton(σ([1.0, 2.0], E)), lab="d1")
plot!(Singleton(σ([-1.0, -2.0], E)), lab="-d1")
plot!(Singleton(σ([-1.0, 2.0], E)), lab="d2")
plot!(Singleton(σ([1.0, -2.0], E)), lab="-d2")
```

The following formulas hold:

- $\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell) =

- $\rho_{M\mathcal{X}}(\ell) =  for any
Here $M\mathcal{X}$ represents the linear map for $M \in \mathbb{R}^{m \times n}$.
- $\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell) =

The following table summarizes the number of
| Operation                                     | Simplification rule | Cost  |
|-----------------------------------------------|-------|------|
| $\rho_{\mathcal{X} \oplus \mathcal{Y}}(\ell)$ |$\rho_{\mathcal{X}}(\ell) + \rho_{\mathcal{Y}}(\ell)$  |  $1$     |
| $\rho_{M\mathcal{X}}(\ell)$                   |$rho_{\mathcal{X}}(M^T\ell)$|   $2mn$   |
| $\rho_{CH(\mathcal{X}, \mathcal{Y})}(\ell)$   |$\max\{\rho_{\mathcal{X}(\ell)}, \rho_{\mathcal{Y}(\ell)}\}$|  $1$    |

which touches and contains $\mathcal{X}$ . If $\ell$ is of unit length, then
$\rho_{\mathcal{X}}(\ell)$ is the signed distance of $\mathcal{H}_{\ell}$ to the origin.
Evaluating the support function for a set of directions $L ⊆ \mathbb{R}^n$ provides an overapproximation
```math
    \lceil \mathcal{X} \rceil _L = \bigcap_{\ell \in L} \{ x \in \mathbb{R}^n | \ell^T x \leq \rho_{\mathcal{X}}(\ell) \}
```

i.e., $\mathcal{X} ⊆ \lceil \mathcal{X} \rceil _L$. If $L = \mathbb{R}^n$, then $\mathcal{X} = \lceil \mathcal{X} \rceil _L$, so the support function represents any convex set $\mathcal{X}$ exactly. If $L$ is a finite set of directions $L = {\ell_1, . . . , \ell_m}$, then $\lceil \mathcal{X} \rceil _L$ is a polyhedron.


## Taylor models

```@example
#=
using ReachabilityAnalysis, Plots

f(x) = -6x^3 + (13/3)x^2 + (31/3)x
dom = -3.5 .. 3.5

plot(f, -3.5, 3.5, lab="f", xlab="x")

x = Taylor1(5)
set_taylor1_varname("x")
f(x)

rem = 0 .. 0
x0 = 0.0
dom = -3.5 .. 3.5
tm = TaylorModel1(f(x), rem, x0, dom)
=#
```

## Approximations

Given a set ``\mathcal{X} \subseteq \mathbb{R}^n``, any subset of ``\mathcal{X}``
is said to be an *underapproximation*. Conversely, any set containing ``\mathcal{X}``
is said to be an *overapproximation*. In this section we recall the definitions
and give some examples of the basic set operations commonly used to construct
reachability algorithms. Such operations are not only required to propagate
reachable sets

### Conversion

### Overapproximation

### Underapproximation

### Hausdorff distance

We begin with the notion of Hausdorff distance which, as a way to *measure*
(in the informal sense) the distance between sets, constitutes a practical tool to
quantify the quality of an approximation.

```math
  d_H(\mathcal{X}, \mathcal{Y}) = \max \left( \sup_{x \in \mathcal{X}}\inf_{y \in \mathcal{Y}} \Vert x - y \Vert, \sup_{y \in \mathcal{Y}}\inf_{x \in \mathcal{X}} \Vert x - y \Vert \right)
```



## References
