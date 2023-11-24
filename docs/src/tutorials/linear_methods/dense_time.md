```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

## Set propagation in dense time

In this section we transition from the discrete time problem considered before
to a problem in which time is assumed to be a continuously varying quantity.
Thus, our next "toy" problem is the system of ordinary differential equations
that describe a simple, undamped, harmonic oscillator:

```math
\left\{ \begin{aligned}
x'(t) &= y(t) \\
y'(t) &= -x(t)
\end{aligned} \right. \qquad \textrm{subject to } (x(0), y(0)) \in X_0 = [0.8, 1.2] \times [0.8, 1.2] \subset \mathbb{R}^2
```
We are interested in solving this problem for all $t \in [0, T]$ with $T = 2\pi$.

The problem studied in the previous section is the "discrete analogue" of the purely continuous system of differential equations defined above.

The matrix $M(\omega t)$ describes the solution of a mathematical model called *simple harmonic oscillator* with natural frequency $\omega = 1$, and the **analytic solution** is

```math
\left\{ \begin{aligned}
x(t) &= x_0 \cos~ t + y_0 \sin~ t \\
y(t) &= -x_0 \sin t + y_0 \cos t
\end{aligned} \right.
```

Observe that $(x(t), y(t))^T$ is just the matrix $M(t)$ applied to the initial state $(x_0, y_0)^T$.


Define the invariant $G: \{(x, y) \in \mathbb{R}^2: x ≥ 1.3 \}$.


```@example dense_propagation
using ReachabilityAnalysis
using ReachabilityAnalysis: center

import Plots: plot, plot!, xlims!, ylims! # hide

import Random # hide
Random.seed!(1117) # hide

# initial states
X0 = BallInf(ones(2), 0.2)

# rotation matrix
M(θ) = [cos(θ) sin(θ); -sin(θ) cos(θ)]

# analytic solution to visualize the solution at intermediate times
function analytic_sol(u0, tt)  # hide
    sol = [M(t) * u0 for t in tt]  # hide
    xcoords = getindex.(sol, 1)  # hide
    ycoords = getindex.(sol, 2)  # hide
    return xcoords, ycoords  # hide
end  # hide

samples = sample(X0, 200, include_vertices=true)
tt = range(0, 2pi/20, length=100)

F = [ReachSet(linear_map(M(ti), X0), ti) for ti in range(0, 2pi, step=2pi/20)]

fig = plot(F[1], vars=(1, 2), lab="X0", c=:blue)
plot!(fig, F[2], vars=(1, 2), lab="Φ * X0", c=:orange)
[plot!(fig, xcoords, ycoords, seriestype=:path, lab="", c=:magenta, ratio=1.) for (xcoords, ycoords) in analytic_sol.(samples, Ref(tt))]
fig

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

### Conservative time discretization



```@example dense_propagation
fig = plot()

plot!(fig, convexify(F[1:2]), vars=(1, 2), ls=:dash, lw=3.0, c=:green, lab="CH(X0, Φ*X0)")
plot!(fig, F[1], vars=(1, 2), lab="X0", c=:blue)
plot!(fig, F[2], vars=(1, 2), lab="Φ*X0", c=:orange)
[plot!(fig, xcoords, ycoords, seriestype=:path, lab="", c=:magenta, ratio=1.) for (xcoords, ycoords) in analytic_sol.(samples, Ref(tt))]

fig

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

```@example dense_propagation
xlims!(1.0, 1.5)
ylims!(0.8, 1.2)
```

```@raw html
<div style="background: ghostwhite;
            font-size: 16px;
            padding: 10px;
            border: 1px solid lightgray;
            margin: 10px;
            font-family: helvetica;">

<b>Proposition 1.</b> Conservative time discretization for linear homogeneous systems [[LGG09]](@ref).

Let $x' = Ax$ with initial states $x(0) \in X_0$. Then the convex set

$$
\Omega_0 := CH(X_0, \Phi X_0 \oplus B_r),
$$

where $\Phi = e^{A \delta}$ is the state transition matrix given a step size
$\delta > 0$ and $B_r$ is the infinity norm ball centered at the origin and
radius $r = (e^{\delta \Vert A \Vert } - 1 - \delta \Vert A \Vert)\Vert X_0\Vert$,
satisfies $R^e([0, \delta], X_0) \subseteq \Omega_0$, where
$R^e([0, \delta], X_0)$ are the exact reachable states between time $0$ and $\delta$.

Moreover, $\Omega_0$ converges to the true reachable states when the step size
$\delta \to 0$ (in Hausdorff distance).

</div>
```

Algorithms implementing conservative time discretization can be used from the
[`discretize(ivp::IVP, δ, alg::AbstractApproximationModel)`](@ref) function.
Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.
This function receives three inputs: the initial value problem (``ivp`) for a
linear ODE in canonical form, (e.g. the system returned by `normalize`);
the step-size (`δ`), and the algorithm (`alg`) used to compute the approximation model.
Do `subtypes(ReachabilityAnalysis.AbstractApproximationModel)` to see the
available approximation models.

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
