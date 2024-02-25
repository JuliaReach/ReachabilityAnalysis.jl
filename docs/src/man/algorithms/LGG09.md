```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Support-function based method (LGG09)

## Method

In the following subsections we outline the method of [[LGG09]](@ref) to solve
linear set-based recurrences using support functions, first the homogeneous case
and then the inhomogeneous case without wrapping effect.
We also discuss the special case of real eigenvalues.

## Homogeneous case

Consider the set-based recurrence

```math
X_{k+1} = Φ X_k,\qquad 0 ≤ k ≤ N
```
where ``Φ ∈ \mathbb{R}^{n\times n}`` and ``X_0 ⊆ \mathbb{R}^n`` are given.
By unrwapping the recurrence, ``X_k = Φ^k X_0`` for all ``k ≥ 0``. Let ``d ∈ \mathbb{R}^n`` be a given
*template direction*. Using the property of support functions ``ρ(d, A X) = ρ(A^T d, X)``
for any matrix ``A`` and set ``X``, we have that

```math
ρ(d, X_k) = ρ(d, Φ^k X_0) = ρ((Φ^T)^k d, X_0).
```
In this way we are able to reason with the sequence ``\{X_0, X_1, X_2, …, X_N\}``
by evaluating the support function of the initial set ``X_0`` along the directions
``\{d, Φ^T d, (Φ^T)^2 d, …, (Φ^T)^N d\}``.


## Inhomogeneous case

The inhomogeneous case generalizes the previous case by taking, at each step,
the Minkowski sum with an element from the sequence ``\{V_0, V_1, V_2, …, V_N\}``:

```math
X_{k+1} = Φ X_k \oplus V_k,\qquad 0 ≤ k ≤ N.
```
Let us write such recurrence in the unrapped form,

```math
\begin{aligned}
\quad X_1 &= Φ X_0 \oplus V_0 \\[1mm]
\quad X_2 &= Φ X_1 \oplus V_1 = Φ^2 X_0 \oplus Φ V_0 \oplus V_1 \\[1mm]
\quad X_3 &= Φ X_2 \oplus V_2 = Φ^3 X_0 \oplus Φ^2 V_0 \oplus Φ V_1  \oplus V_2 \\[1mm]
\quad &\vdots \\[1mm]
\quad X_k &= Φ^k X_0 \oplus \left( \bigoplus_{i=0}^{k-1} Φ^{k-i-1} V_i \right)
\end{aligned}
```
where the big Minkowski sum is just an abbreviation for
``Φ^{k-1} V_0 \oplus Φ^{k-2} V_1 \oplus Φ^{k-3} V_2 \oplus … \oplus Φ V_{k-2} \oplus V_{k-1}``.

Let ``d ∈ \mathbb{R}^n`` be a given template direction. Using the additive property of
support functions, ``ρ(d, X \oplus Y) = ρ(d, X) + ρ(d, Y)`` for any sets ``X`` and ``Y``,
we have that

```math
\begin{aligned}
\quad ρ(d, X_1) &= ρ(Φ^T d, X_0) + ρ(d, V_0) \\[1mm]
\quad ρ(d, X_2) &= ρ((Φ^T)^2 d, X_0) + ρ(Φ^T d, V_0) + ρ(d, V_1) \\[1mm]
\quad ρ(d, X_3) &= ρ((Φ^T)^3 d, X_0) + ρ((Φ^T)^2 d, V_0) + ρ(Φ^T d, V_1) + ρ(d, V_2) \\[1mm]
\quad &\vdots \\[1mm]
\quad ρ(d, X_k) &= ρ((Φ^T)^k d, X_0) + \sum_{i=0}^{k-1} ρ( (Φ^T)^{k-i-1} d,  V_i).
\end{aligned}
```
In a similar fashion to the homogeneous case, the method allows to efficiently reason
about the the sequence ``\{X_0, X_1, X_2, …, X_N\}`` by evaluating the support
function of the initial set ``X_0`` and the input sets ``\{V_k\}_k`` along the directions
``\{d, Φ^T d, (Φ^T)^2 d, …, (Φ^T)^N d\}``. Implementation-wise, we update
two sequences, one that accounts for the homogeneous term, and another
sequence that accounts for the effect of the accumulated inputs.

## Implementation details

The reach-set representation used is a [`TemplateReachSet`](@ref), which stores the
directions used (vector of vectors) and the support function evaluated at each direction
(matrix, see below). The set representation, `set(R::TemplateReachSet)`, is either a polyhedron in constraint form
(`HPolyhedron`), or a polytope (`HPolytope`) if the directions are bounding, i.e.
the template directions define a bounded set.

The computed support function values can accessed directly through the field
`sf::SN` of each template reach-set. Here `sf` is an array view of type `::Matrix{N}(undef, length(dirs), NSTEPS)`:
each row corresponds to one of the template directions and each column corresponds to a fixed iteration index ``k ≥ 0``.

If you use directions from the canonical basis of ``\mathbb{R}^n``, it is recommended to define `LazySets.Arrays.SingleEntryVector`
or "one-hot" arrays as they are commonly called, because there are methods that dispatch on such type of arrays efficiently.

## Parallel formulation

The support functions of the sequence ``\{X_k\}_k`` along different directions can be
computed in parallel. Indeed, if ``d_1`` and ``d_2`` are two given template directions, two different processes
can independently compute ``ρ(d_1, X_k)`` and ``ρ(d_2, X_k)`` for all ``k = 0, 1, …, N``
using the methods described above. Once both computations have finished, we can store
the resulting support functions in the same array. Use the flag `threaded=true` to
use this method.

Implementation-wise the function `_reach_homog_dir_LGG09!` spawns different threads
which populate the matrix `ρℓ::Matrix{N}(undef, length(dirs), NSTEPS)` with the computed
values. Hence each thread computes a subset of distinct rows of `ρℓ`.

## Real eigenvalues

If the spectrum of the state matrix only has real eigenvalues, the sequence of
support functions can be computed efficiently if we work with a template
consisting of eigenvectors of ``Φ^T``. This idea is described in [[LGG09]](@ref)
and we recall it here for convenience.

The method stems from the fact that if ``(λ, d)`` is an eigenvalue-eigenvector
pair of the matrix ``Φ^T``, with ``λ ∈ \mathbb{R}``, then
``Φ^T d = λ d``, and if we apply ``Φ^T`` on both sides of this identity, we get
``(Φ^T)^2 d = Φ^T (Φ^T d) = Φ^T(λ d) = λ^2 d``.
In more generality, it holds that ``(Φ^T)^k d  = λ^k d`` for all ``k ≥ 1``.
Applying this relation to the support function recurrence described above, we get
for the general inhomogeneous and possibly time-varying inputs case:

```math
ρ(d, X_k) = ρ(λ^k d, X_0) + \sum_{i=0}^{k-1} ρ(λ^{k-i-1} d, V_i).
```
To further simplify this formula, we analyze different cases of ``λ``.
If ``λ = 0``, then ``ρ(d, X_k) = ρ(d, V_k)`` for all ``k ≥ 1``, so we focus
on either ``λ`` being positive or negative. To further simplify the computation
of ``ρ(d, X_k)``, we can use the property ``ρ(λ d, X) = λ ρ(d, X)``
if ``λ ≥ 0``. We now consider the cases ``λ > 0`` and ``λ < 0``.

**Case ``λ > 0``.** Then ``λ^k > 0`` for all ``k ≥ 1``,
and

```math
ρ(d, X_k) = λ^k ρ(d, X_0) +  \sum_{i=0}^{k-1} λ^{k-i-1} ρ(d, V_i).
```
We are left with evaluating the support function only at ``ρ(d, X_0)``  and ``ρ(d, V_i)``
to construct the full sequence ``\{ρ(d, X_k)\}_{k}``. Moreover, if the ``V_i``'s are constant
we can extract them from the right-hand side sum and use that
```math
\sum_{i=0}^{k-1} λ^{k-i-1} = 1 + λ + … + λ^{k-1} = \dfrac{1 - λ^k}{1 - λ}.
```

**Case ``λ < 0``.** Since ``λ^k = (-1)^k (-λ)^k`` and ``λ < 0``, then
``λ^k`` is positive if ``k`` is even, otherwise it is negative. So we can write:

```math
ρ(d, X_k) = (-λ)^k ρ((-1)^k d, X_0) + \sum_{i=0}^{k-1} (-λ)^{k-i-1} ρ((-1)^{k-i-1} d, V_i).
```
The main difference between this case and the previous one is that now we have to evaluate
support functions ``ρ(\pm d, X_0)`` and ``ρ(\pm d, V_i)``. Again, simplification takes place
if the ``V_i``'s are constant and such special case is considered in the implementation.
