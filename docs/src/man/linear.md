```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Linear ordinary differential equations


## Conservative time discretization

Algorithms implementing conservative time discretization can be used from the
[`discretize(ivp::IVP, δ, alg::AbstractApproximationModel)`](@ref) function.
Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.
This function receives three inputs: the initial value problem (``ivp`) for a linear ODE in canonical form,
(e.g. the system returned by `normalize`); the step-size (`δ`), and the algorithm
(`alg`) used to compute the approximation model.


```@example
subtypes(ReachabilityAnalysis.AbstractApproximationModel)
```

The output of a discretization is a new initial value problem of a discrete system.
Different approximation algorithms and their respective options are described
in the docstring of each method, e.g. [`Forward`](@ref).


Initial-value problems considered in this function are of the form

```math
x' = Ax(t) + u(t),\\qquad x(0) ∈ \\mathcal{X}_0,\\qquad (1)
```
and where ``u(t) ∈ U(k)`` add where ``\\{U(k)\\}_k`` is a sequence of sets of
non-deterministic inputs and ``\\mathcal{X}_0`` is the set of initial
states. Recall that this initial-value problem is called homogeneous whenever `U`
is the empty set. Other problems, e.g. ``x' = Ax(t) + Bu(t)`` can be brought
to the canonical form with the function [`normalize`](@ref).

The initial value problem returned by this function consists of a set discretized
(also called *bloated*) initial states ``Ω₀``, together with the coefficient matrix
``Φ = e^{Aδ}`` and a transformed sequence of inputs if ``U`` is non-empty.

Two main variations of this algorithm are considered: dense time case and
discrete time case.

- In the dense time case, the transformation is such that the trajectories
of the given continuous system are included in the computed flowpipe of the
discretized system. More precisely, given a step size ``δ`` and the system (1)
conservative set-based discretization function computes a set, ``Ω₀``, that
guarantees to contain all the trajectories of (1) starting at any ``x(0) ∈ \\mathcal{X}_0``
and for any input function that satisfies ``u(t) ∈ U(1)``, for any ``t ∈ [0, δ]``.
If ``U`` is time-varying, this function also discretizes the inputs for ``k ≥ 0``.

- In the discrete time case, there is no bloating of the initial states and the
input is assumed to remain constant between sampled times. Use the algorithm
`NoBloating()` for this setting. If ``U`` is time-varying, this function also discretizes
the inputs for ``k ≥ 0``.

There are algorithms to obatin such transformations, called *approximation models*
in the technical literature. For references to the original papers, see the
docstring of each concrete subtype of `AbstractApproximationModel`.
