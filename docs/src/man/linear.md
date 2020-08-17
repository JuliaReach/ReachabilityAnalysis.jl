# Linear ordinary differential equations

## An illustrative example

The textbook example of a damped spring-mass system corresponds to the following linear ODE:

```math
\dfrac{d^2 x}{dt^2} + 2 \xi \omega_0\dfrac{d x}{dt} + \omega_0^2 x = 0
```

Then we consider the [spring-mass system](https://en.wikipedia.org/wiki/Simple_harmonic_motion#Examples),
a second order linear ODE studied in introductory physics courses. For that system
we show an example of how to compute and project the flowpipe, and then plot the
variables of interest.

The solution process mainly consists of three steps:

(1) **Formulating** the mathematical problem, in the form of an initial-value problem (IVP)
    with possibly uncertain initial states or inputs.

(2) **Solving** the IVP, either with the default algorithm or specifying the algorithm
    and some of its options.

(3) **Extracting** the results, either to visualize with a plot, or to project onto
    the relevant variables for further study.

Below we give some further details of this solution step for the simple scalar equation
presented above.



## Set-based recurrences

Set-based methods can be applied for verifying safety properties of dynamical systems.
In this section we overview the set-based reachability problem for affine dynamical systems,

```math
x'(t) = Ax(t) + Bu(t), \qquad (1)
```
where ``x(t) \in \mathbb{R}^n`` is the state-space vector and ``A`` and ``B`` are
matrices of the appropriate dimensions. We consider an initial state that can be
any point in a given set ``x(0) \in \mathcal{X}_0`` and ``u(t) \in \mathcal{U} \subseteq \mathbb{R}^{m}``
is a non-deterministic input. Both the initial set and the set of input functions
are assumed to be compact and convex. Here "non-deterministic input" refers to any bounded, measurable,
input functions $u(t)$ which are included in the given set ``\mathcal{U}``.
An extension to time-varying input sets is described later on.

It is usual that one is interested in computing *observable* outputs,

$$
y(t) = Cx(t) + Du(t),\qquad (2)
$$
where $C$ and $D$ are matrices of the appropriate dimension. In mathematical systems theory, equations (1-2) are known as *linear time-invariant* or just LTI systems.

!!! note "Handling non-convex sets"
    There is increasing support for non-convex set represenations within `Lazyset.jl`.
    One approach to handle non-convex initial sets is to cover such sets with one or more convex sets;
    if several sets are used, the split initial sets API can be used and the different flowpipes are computed
    in parallel.

In this context, set- based recurrence relatinos of the form

$$
\mathcal{X}(k+1) = \Phi\mathcal{X}(k) \oplus \mathcal{V}(k),\qquad k = 0,1 ,\ldots, N,\qquad (3)
$$
arise naturally.

if you are familiar with ODEs . . . Euler discretization...

## Approximation model
