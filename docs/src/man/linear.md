```@meta
DocTestSetup  = quote
    using MyPackage
end
CurrentModule = ReachabilityAnalysis
```

# Linear ordinary differential equations

## Introduction

Consider a ball thrown upwards from the ground . . .

The one-dimensional scalar ODE with exponential decay,

```math
x'(t) = -x(t),\qquad t ∈ [0, T]
```
is among the first studied in elementary calculus courses. We can compute the
flowpipe as follows:

```@example linear_scalar
using ReachabilityAnalysis, Plots

# define an initial-value problem
prob = @ivp(x' = -x, x(0) ∈ 0.45 .. 0.55)

# solve it
sol = solve(prob, T=4.0)

# plot the solution, where the index 0 corresponds to the "time" variable
fig = plot(sol, vars=(0, 1), label="Flowpipe", xlab="t", ylab="x(t)")

# plot some trajectories
x0vals = [0.45, 0.50, 0.55]
[plot!(fig, t -> x0*exp(-t), xlims=(0, 4), label="Analytic with x(0) = $x0") for x0 in x0vals]
fig
```
The figure shows the flowpipe of the initial-value problem
$x'(t) = -x(t)$, $0 ≤ t ≤ T = 4.0$, with an initial condition
in the set $x(0) ∈ [0.45, 0.55]$. Each blue box is called a *reach-set*,
as it represents the set of states the system can reach for a given time-span.
The union of all reach-sets is called the *flowpipe*. The time-span of the first
three reach-sets is

```@example linear_scalar
[tspan(R) for R in sol[1:3]]
```
Note that the time spans are equally separated, because the algorithm used
has a constant time-step.

-- obtained with a default algorithm choice --

In this simple case we know that for an initial point $x_0 \in \mathbb{R}$,
the solution is $x(t) = x_0 e^{-t}$, so we plotted some trajectories
as well.

!!! note "Performance tip"
    It is often instructive to compare runtimes between different solution methods,
    even though reachability vs. traditional numerical integration is truly
    comparing [apples and oranges](https://en.wikipedia.org/wiki/Apples_and_oranges):
    in the former case, there is a notion of an infinity of trajectories being computed,
    as the object of study (the exact flowpipe) is naturally represented with sets;
    in the latter case, no numerical overapproximation is required, etc.
    That said, the computation time for the exponential decay problem with step sizes
    `1e-1`, `1e-2` and `1e-3` is `12.504 μs`, `33.198 μs` and `232.155 μs` respectively,
    while the time taken to evaluate the analytic solution on the considered time
    span is `558.065 ns`, `4.322 μs` and `39.907 μs` respectively. These experiments were ran on Julia v1.5 and a standard laptop (Intel i7 CPU@3.10GHZ).

Observing the region covered by the flowpipe we see that it strictly includes
that covered by the trajectories at the extreme values,
$x(0) = 0.45$ and $x(0) = 0.55$, but it also includes some additional area:
there are small triangular regions of the flowpipe which *do not correspond to any
trajectory* of the system. This is just the result of computing reach-set
overapproximations: we know for sure that *all trajectories* are included, though
the set computed may also contain other regions of state-space which are spurious.
To overcome this effect there are two common alternatives:

- Either reduce the time-step of the algorithm, e.g. try passig `alg=INT(δ=1e-2)`
  to the `solve` function to wit that the flowpipe is much tighter, or
- work with under-approximation, which basically suffer from the dual problem of
  possibly missing some trajectories, but not adding extra ones.

It is worth noting that theoretical estimates on the [Hausdorff distance](@ref) between
the exact flowpipe and the computed flowpipe for a given time-step are known, but
are often not very useful in practice since they tend to provide too coarse bounds.

The specification of constraints on the evolution can be given in terms of
*state invariants.* For example, suppose that we restrict to the set $X ≥ 0.1$.

```@example linear_scalar
f(δ) = solve(prob, T=4.0, alg=INT(δ=δ))

# plot over the full time span
fig = plot(xlab="t", ylab="x(t)")
plot!(fig, sol, vars=(0, 1), label="", xlab="t", ylab="x(t)" , lw=0.1, c=:orange)
plot!(fig, HalfSpace([0.0, -1.0], -0.1), alpha=.2, c=:green, lab="X")
plot!(fig, f(1e-2), vars=(0, 1), label="", , lw=0.0, c=:blue)
x0sample = [0.45, 0.50, 0.55]
[plot!(fig, t -> x0*exp(-t), xlims=(0, 4), label="") for x0 in x0sample]
plot!(fig, Hyperrectangle(low=[1.09, 0.05], high=[2.0, 0.15]), alpha=0.2, lw=2.0, c=:black, linestyle=:dot)

# zoomed plot
figz =


```

For illustration, suppose that we are interested in the behavior of trajectories
escaping the invariant. Define the *guard set* as the (one-dimensinal) hyperplane
$\partial X := \{x \in \mathbb{R} : x = 1\}$, i.e. it is in this case the border
of our invariant $X$. We can analytically compute the time at which the trajectory
with initial condition $x(0) = x_0$ exits the invariant
by solving for $t$ in the analytic equation, $x(t) = x_0 e^{-t}$, giving
$t = -log_e(0.1/x_0)$. Evaluating for the extremal values of $$\mathcal{X}_0 = [0.45, 0.55]$$,
we get:

```@example linear_scalar
xref = -log(0.1 / 0.45) .. -log(0.1 / 0.55)
```
For comparison, let's compare the *reference interval* `xref` computed above with the
bracketing interval obtained with reachability analysis, for different values of the
step size ``δ``.

```@example linear_scalar
using ReachabilityAnalysis: relative_error

using PrettyTables, BenchmarkTools

steps = 10 .^ range(-2, -7, length=5)
data = Matrix{Any}(undef, length(steps), 4)

for (i, δ) in enumerate(steps)
    bench = @benchmark sol($δ)
    btime = minimum(bench.times) * 1e-6 # ns to ms

    local sol = f(δ)
    idx = findall(R -> [0.1] ∈ R, sol)
    x = tspan(sol[idx])
    rerr = relative_error(x, xref)

    data[i, :] = Any[δ, x, rerr, btime]
end

tags = ["δ", "Bracketing interval", "Relative error (%)", "Runtime (ms)"]
pretty_table(data, tags; formatters = ft_printf("%5.4f", 4))
```

!!! note "Performance tip"
    In this simple example we used a fixed step-size $\delta$. It is
    an interesting exercise to write an algorithm that dynamically chooses the step-size,
    by refining only the subset of the flowpipe that reaches the boundary. On one hand,
    such dynamic algorithm shall converge faster; on the other hand, it has the increase
    in complexity of having to handle several flowpipes -- an exponentially
    increasing number of flowpipes -- hence, taking the convex hull (called
    *convexification* in this library) shall be considered.

The examples studied so far illustrate that the solution process mainly consists of three steps:

(1) **Formulating** the mathematical problem, in the form of an initial-value problem
    with possibly uncertain initial states, inputs or both. This is usually done with the
    help of the `@ivp` macro, that recognises the type of system and creates an appropriate
    Julia type. In later sections it is shown how to formulate problems
    in which the parameters of the system's dynamics are also uncertain.

(2) **Solving** the problem, either with the default algorithm coice or specifying the algorithm
    and some options, through the `alg=...` keyword algorithm to `solve`.
    Typical algorithms are the zonotope based algorithm
    [GLGM06](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/lib/algorithms/GLGM06/)
    and the support function based algorithm
    [LGG09](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/lib/algorithms/LGG09/).
    A more advanced alternative not discussed in this section is the decomposition approach
    [BFFPSV18](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/lib/algorithms/BFFPSV18/).

(3) **Extracting** the results, either to visualize with a plot, or to project onto
    the relevant variables for further study. Usually we load `Plots` and use the command
    `plot(sol, vars=(...))`. Please note that the time variable is associated with `0`,
     i.e. to plot the flowpipe associated to variable ``x_5(t)`` as a function of time
     one would write `vars=(0, 5)`. Plots in phase space, like ``x_1(t)`` vs. ``x_2(t)``
     are done passing ``vars=(1, 2)``.

In the rest of this section we study in more detail each of these steps. First,
we illustrate the method [GLGM06](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/lib/algorithms/GLGM06/).
We also show how to efficiently compute reachable states in high dimension using
[LGG09](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/lib/algorithms/LGG09/).
Finally, we consider in more detail the class of linear dynamics equations of motion,
which typically arise in the spatial discretization of partial differential equations (PDEs).

## Set-based recurrence

In this section of the manual we focus on initial-value problems for linear differential
equations of the form

```math
x'(t) = Ax(t) + u(t), \qquad x(0) \in \mathcal{X}_0, x(t) ∈ \mathcal{X}, u(t) \in \mathcal{U}(t),\qquad (1)
```
where ``x(t) \in \mathbb{R}^n`` is the state vector, ``A \in \mathbb{R}^{n \times n}``
is the state (or coefficients) matrix, ``\mathcal{X}_0`` is the set of initial states,
``\mathcal{U}(t) \subseteq \mathcal{R}^m`` is the (possibly time-varying) set of uncertain inputs, and
``\mathcal{X} \in \subseteq \mathcal{R}^n`` is the state invariant.
We also consider observable outputs,

```math
  y(t) = Cx(t) + Du(t),\qquad (2)
```
where $C$ and $D$ are matrices of appropriate dimension. In mathematical systems parlance, Eqs. (1) and (2) above encompasses the class of state-constrained linear time-invariant (LTI) systems with uncertain initial states and uncertain input sets. With respect to the time domain, we consider a finite time span, typically of the form ``[0, T]``, where ``T`` is called the time horizon. Applications that require reasoning about an infinite (unbounded, eg. ``[0, \infty]``) time horizon can make use of so-called invariant set computations; see e.g. the [Van der Pol](@ref) model for a practical example.

!!! tip "Extension tip"
    Please note that systems of the form ``x'(t) = Ax(t) + Bu(t)`` can be brought
    to the canonical form ``x'(t) = Ax(t) + v(t)`` by constraining ``v(t)`` to be
    in the set ``B \mathcal{U}``. Such transformation is implemented in the function
    [normalize](@ref).

Although many systems can be formulated in this manner, arguably most systems of interest
in science and engineering have additional complexities, e.g. the state matrix `A`
can depend on time, the inputs can be state-dependent, or more generally right-hand
side can contain powers or other nonlinear functions of the state `x`, etc. Moreover,
in several situations of interest the coefficients of ``A`` are only known approximately.
On the other hand, studying linear systems is of central interest for our purposes

For example, . . .

Although many differential equations can be formulated in such form,

At its core, any reachability method in the class considered aims at solving the
set-based recurrence

```math
X_{k+1} = \Phi X_k \oplus V_k,\qquad k = 1, 2, …, N
```
In order to solve this recurrence efficiently and accurately, one has to make an appropriate
choice of set representation, i.e. the way in which each ``X_k`` is stored (as well as ``V_k``),
and this decision has an impact of how well the linear map and the Minkowski, both
set operations involved above, are computed. In addition to the set representation,
one has to choose appropriate overapproximation functions, because the complexity
of each operation typically increases with ``k``, unless one simplifies or overapproximates
the set with another set of the same family but with simpler complexity.

In the rest of this section we show how to solve initial value problems with sets of initial
conditions. First we consider the simple *scalar* equation $x'(t) = -x(t)$,
where $x(0)$ may be any point in a given initial interval $\mathcal{X}_0$.
Then we consider the [spring-mass system](https://en.wikipedia.org/wiki/Simple_harmonic_motion#Examples),
a second order linear ODE studied in introductory physics courses. For that system
we show an example of how to compute and project the flowpipe, and then plot the
variables of interest.

Some more advanced topics than those presented here can be found in the section
[Exploiting structure](@ref).

## Zonotope methods

```@example linear_methods
using ReachabilityAnalysis, Plots, Colors

red, green, blue, purple = Colors.JULIA_LOGO_COLORS

A = [0 1; -1 0.]
B = Ball2([0.0, 1.0], 0.1)
prob = @ivp(x' = A * x, x(0) ∈ B)
ω = 2π

sol = solve(prob, tspan=(0.0, 0.9ω), alg=GLGM06(δ=0.2, max_order=10));

plot(sol, vars=(1, 2), lw=0.0, ratio=1.)
plot!(sol(0.0), vars=(1, 2), c=red, alpha=1.)
plot!(sol(.2ω), vars=(1, 2), c=green, alpha=1.)
plot!(sol(.4ω), vars=(1, 2), c=blue, alpha=1.)
plot!(sol(.75ω), vars=(1, 2), c=purple, alpha=1.)
```

## Support function methods

```@example linear_methods
A = [0 1; -1 0.]
B = Ball2([0.0, 1.0], 0.1)
prob = @ivp(x' = A * x, x(0) ∈ B)
ω = 2π

dirs = PolarDirections(100)
sol = solve(prob, tspan=(0.0, 0.9ω), alg=LGG09(δ=0.2, template=dirs));

plot(sol, vars=(1, 2), lw=0.0, ratio=1.)
plot!(sol(0.0), vars=(1, 2), c=red, alpha=1.)
plot!(sol(.2ω), vars=(1, 2), c=green, alpha=1.)
plot!(sol(.4ω), vars=(1, 2), c=blue, alpha=1.)
plot!(sol(.75ω), vars=(1, 2), c=purple, alpha=1.)
```

```@example linear_methods
plot(B, ratio=1., lab="B")
plot!(sol[1], vars=(1, 2), lab="R[1]")
plot!(sol[2], vars=(1, 2), lab="R[2]")
```

## Working with solutions

Extracting relevant information from flowpipes is the next step after the solution
of the initial-value problem has been found.

For further examples, explore the [Examples](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/man/examples_overview/)
section in the documentation.

 we need to process sets of states
(i) efficiently, and (ii) rigorously. About (i), since linear reachability is often
used as a primitive for other methods (eg. used iteratively for nonlinear reachability,
or to perform different kinds of parameter synthesis), efficiency is one of the key
aspects for scalability. Moreover, high dimensional (eg. problems arising from
discretization of PDEs may have many degrees freedom), to the method requires
to work efficiently in high dimensions as well. If the time scales are small,
flowpipes can contain thousands or millions of reach-sets.


A common requirement in set-based computations is being able to extract relevant

Once the flowpipe has been computed, we are interested in extracting information
from it.

!!! note "Performance tip"
    It is possible to check safety properties "on the loop", i.e. without actually
    computing the whole flowpipe. The key idea is to compute the flowpipe lazily,
    and only take concrete operations for the combination of variables that enter in
    the property that we want to check.

To work further we consider a less trivial example in this section.

A 28-dimensional controlled helicopter model is now available here. It illustrates the precision and speed of SpaceEx. A SpaceEx flowpipe with box constraints in this 28-dimensional model takes about 15 seconds, and corresponds to 228 corner case simulations -- that's more than 268 million simulation runs.


In the following example we consider a spring-mass system which is a linear ODE
with two-degrees of freedom, ``x''(t) = -x(t)``, or:

```math
x'(t) = v(t),\qquad v'(t) = -x(t), \qquad t \in [0, T]
```

## Linear dynamics equations of motion

In this section we consider systems of linear differential equations of second order,

```math
    Mx''(t) + Cx'(t) + Kx(t) = f(t),\qquad t \in [0, T], \qquad (1)
```
where ``M``, ``C`` and ``K`` are the mass, viscosity (or damping) and stiffness matrices respectively, and ``f(t)`` is the vector of externally applied loads. A noteworthy application of such systems is the equation of equilibrium governing the linear dynamic response of a system of finite elements (FEM). Here, `x(t)`, `v(t) = x'(t)` and `a(t) = x''(t)` are the displacement, velocity, and acceleration vectors of the finite element assemblage. Moreover, due to physical considerations it is often the case that ``M``, ``C`` and ``K`` are symmetric; ``M`` is positive-definite and ``C`` and ``K`` are positive-semidefinite. We refer to [^BAT06] for further details on the FEM context.

The initial-value problem for Eq.~(1) is consists of finding a displacement ``x(t)`` satisfying Eq.~(1) and given initial data,

```math
x(0) \in X_0,\qquad v(0) ∈ V_0.
```

In the rest of this section we formulate and solve problem (1) in dense time for uncertain initial conditions and input forces using linear reachability methods. We assume that ``M`` is invertible, in which case Eq.~(1) can be transformed into a system of first order ODEs by multiplying with ``M^{-1}`` on the left and introducing the auxiliary vector ``[x(t),~v(t)]``.

### Free oscillations

Our first example is a damped oscillating system without forcing term,
```math
    x''(t) + 0.5~x'(t) + 4x(t) = 0, \qquad x(0) ∈ [0.7 .. 1.3], v(0) ∈ [0.2 .. 0.8].
```
The implementation and solution is straightforward.

```@example second_order_damped
using ReachabilityAnalysis, Plots

# x'' + 0.5x' + 4x = 0
sys = SecondOrderLinearContinuousSystem(hcat([1.0]), hcat([0.5]), hcat([4.0]))

X0 = Interval(0.7, 1.3)
V0 = Interval(0.2, 0.8)
prob = @ivp(sys, x(0) ∈ (X0, V0))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01));

plot(sol, vars=(0, 1), lw=.1, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), lw=0.1, xlab="time", lab="v(t)")
```
We used the zonotope-based algorithm `GLGM06` with step-size`δ=0.01` and
plotted the flowpipe for $x(t)$ and $v(t)$ variables. The time horizon was chosen large enough
as to show the damping effect.

As a simple illustration of working with solutions, suppose that we are interested
in the maximum and minimum values of the velocity ``v(t)`` of this system for
``t ≥ 10s``. We can filter the solution for the time interval between ``10s`` and
``20s`` with the command ``sol(10 .. 20)``. The maximum value of the velocity
is obtained by computing the support function of the flowpipe along direction
``[0.0, 1.0]``; a similar argument works for the minimum of the solution, but
with a ``-`` sign in front, as shown below.

```@example second_order_damped
vmax = ρ([0.0, 1.0], sol(10 .. 20))
vmin = -ρ([0.0, -1.0], sol(10 .. 20))

P = HPolyhedron([HalfSpace([0.0, 1.0], vmax),    # v <= vmax
                 HalfSpace([0.0, -1.0], -vmin),  # v >= vmin
                 HalfSpace([-1.0, 0.0], -10.0)]) # t >= 10

plot!(x -> 0.25, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(x -> -0.25, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(P, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:lightgreen)
```
Where we have plotted the polyhedron with max and min bounds of the flowpipe.
We have thus *proved* that ``-0.25 ≤ v(t) ≤ 0.25`` for all ``t ∈ [10, 20]`` and
for all possible values of ``x(0) ∈ X_0``, ``x'(0) ∈ V_0``.

!!! tip "Extension tip"
    To be mathematically rigorous, we would have to use set types with coefficients
    which are intervals (instead of floating-point numbers) and use the
    `IntervalArithmetic.jl` package to perform absolutely all the intermediate computations.

Let us now analyze the behavior of the solution with the damping coefficient.
If the damping coefficient is higher, then the oscillations decrease more rapidly.
Conversely, by decreasig the damping coefficient we can expect that eventually,
the property ``-0.25 ≤ v(t) ≤ 0.25`` for ``t ∈ [10, 20]`` is no longer satisfied.

```@example second_order_damped
# x'' + cx' + 4x = 0
sys_param(c) = SecondOrderLinearContinuousSystem(hcat([1.0]), hcat([c]), hcat([4.0]))
X0 = Interval(0.7, 1.3)
V0 = Interval(0.2, 0.8)
ivp(c) = @ivp(sys_param(c), x(0) ∈ (X0, V0))
sol_param(c, δ) = solve(ivp(c), tspan=(0.0, 20.0), alg=GLGM06(δ=δ));

cmin, cmax = 0.2, 0.5
cvals = range(cmin, cmax, step=0.01)
sol_tot2(δ) = sol_param.(cvals, δ);
max_v2(δ) = [ρ([0.0, 1.0], fp(10 .. 20)) for fp in sol_tot2(δ)]
min_v2(δ) = [-ρ([0.0, -1.0], fp(10 .. 20)) for fp in sol_tot2(δ)]

fig = plot(xlab="c", ylab="Velocity bounds")

plot!(fig, cvals, max_v2(0.2), lab="v max (δ=0.2)", lw=2.0, c=:green, linestyle=:dashdot)
plot!(fig, cvals, min_v2(0.2), lab="v min (δ=0.2)", lw=2.0, c=:blue, linestyle=:dashdot)

plot!(fig, cvals, max_v2(0.01), lab="v max (δ=0.01)", lw=2.0, c=:green)
plot!(fig, cvals, min_v2(0.01), lab="v min (δ=0.01)", lw=2.0, c=:blue)

plot!(fig, x -> 0.25, xlims=(cmin, cmax), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(fig, x -> -0.25, xlims=(cmin, cmax), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
```

### Constant forcing term

Consider a simple system for which the governing equations are

```math
    Mx''(t) + Kx'(t) = R,
```
where the mass matrix ``M``, stiffness matrix ``K`` and forcing term ``R`` are defined, respectively, as:
```math
M = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}, \qquad M = \begin{pmatrix} 6 & -2 \\ -2 & 4 \end{pmatrix},\qquad R = \begin{pmatrix} 0 \\ 10 \end{pmatrix}
```
This example is taken from Chapter 9 in [^BAT06], where the anaytic solution for
null initial conditions, ``x(0) = x'(0) = 0``, is shown to be

```math
\begin{pmatrix} x(t) \\ x'(t) \end{pmatrix} =
\begin{pmatrix} \frac{1}{\sqrt{3}} & \frac{1}{2}\sqrt{\frac{2}{3}} \\  \frac{1}{\sqrt{3}} & -\sqrt{\frac{2}{3}} \end{pmatrix}\begin{pmatrix} \frac{5}{\sqrt{3}}(1 - \cos(t\sqrt{2}) \\ 2\sqrt{\frac{2}{3}}(-1 + \cos(t\sqrt{5}) \end{pmatrix}
```
Let us first compare the analytic solution with the reachability computation
in discrete time (i.e. without bloating) and singleton initial conditions, obtained
using the `ORBIT` algorithm. This algorithm returns a sequence of singletons that
match the exact solution of the ODE at multiples of the step size. In the same
plot we show the reachability computation with step size ``δ=1-2`` using the
algorithm `GLGM06`.

```@example second_order_damped
using ReachabilityAnalysis, Plots

# Mx'' + Kx = R
M = [2 0; 0 1.]
K = [6 -2; -2 4.]
C = zeros(2, 2)
R = [0, 10.]
sys = SecondOrderAffineContinuousSystem(M, C, K, R)
X0 = zeros(2)
V0 = zeros(2)
prob = @ivp(sys, x(0) ∈ (X0, V0))

# analytic solution
A = [1/√3  (1/2)*√(2/3);
     1/√3      -√(2/3)]
x₁(t) = (5 / √3) * (1 - cos(t*√2))
x₂(t) = (2 * √(2/3)) * (-1 + cos(t*√5))
U(t) = A * [x₁(t), x₂(t)]

fig = plot(xlab="time", ylab="x1(t)")

# solution without bloating
δ = 0.1
sol_orbit = solve(prob, tspan=(0.0, 20.0), alg=ORBIT(δ=0.1))
plot!(fig, sol_orbit, vars=(0, 1), lw=2.0, markershape=:star8)

# solution with bloating
sol_lgg = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(fig, sol_lgg, vars=(0, 1), lw=0.0, lab="GLGM06", c=:green)

# analytic solution
tdom = range(0.0, 20.0, step=0.001)
plot!(fig, tdom, [U(ti)[1] for ti in tdom], lab="Analytic", c=:magenta)

fig
```

The accumulated error is (numerically speaking) zero for the solution in discrete
time.


```@example second_order_damped
# check total error
x = k -> U((k-1) * δ)[1]
x′ = k -> U((k-1) * δ)[2];
```

```@example second_order_damped
sum((set(R).element[1] - x(k))^2 for (k, R) in enumerate(sol_orbit))
```

```@example second_order_damped
sum((set(R).element[2] - x′(k))^2 for (k, R) in enumerate(sol_orbit))
```

Finally, let us also note that by we can solve for a set of initial
conditions, for example, a box around zero, by changing the initial-value problem:

```@example second_order_damped
fig = plot(xlab="time", ylab="x1(t)")
X0box(ε) = BallInf(zeros(4), ε)

# solution with bloating and a box of initial conditions
ε = 1.0
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(fig, sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

ε = 0.5
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

ε = 0.1
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

fig
```

### Nonlinear forcing term

In cases where the forcing term is a nonlinear elementary function,
eg. a combination of trigonometric functions, introducing auxiliary
variables can be used to represent such function.

For instance, if ``f(t) = F \sin (3t)``, let ``u := F\sin(3t)`` and
``g := F \cos(3t)``. Then, it is a simple exercise to see that the system
```math
    x''(t) + 3.5~x'(t) + 4x(t) = F\sin(3t), \qquad x(0) ∈ [0.7 .. 1.3] × [0.7 .. 1.3]
```
is formally equivalent to the following linear system:

```@example second_order_damped

F = 3.0
A = [ 0    1.0     0   0;
     -4    -3.5    1   0;
      0     0      0   F;
      0     0     -1   0]
sys = @system(x' = A * x)

B0 = BallInf(ones(2), 0.3) × Singleton([0.0]) × Interval(0.8, 2.1)
sol = solve(@ivp(sys, x(0) ∈ B0), tspan=(0.0, 15.0), alg=GLGM06(δ=0.01));

plot(sol, vars=(0, 1), c=:magenta, lw=.0, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), c=:green, lw=.0, xlab="time", lab="v(t)")
```

## References

- [^BAT06]: Bathe, Klaus-Jürgen. *Finite element procedures.* Klaus-Jurgen Bathe, 2006.
