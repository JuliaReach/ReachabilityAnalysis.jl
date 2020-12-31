```@meta
DocTestSetup  = quote
    using MyPackage
end
CurrentModule = ReachabilityAnalysis
```

# Linear ordinary differential equations

## Introduction

Consider the following simple physical process: a ball is thrown from the group up.

<<<<<<<<<<<<<


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

## Vertex representation

It is also possible to solve the linear set-based recurrence using the vertex
representation. Such method is implemented in the algorithm [`VREP`](@ref).
In general, using the vertex representation scales worse than other methods such as
zonotopes of support functions. However, sometimes we are either working in low dimensions,
or we want to use initial states of some particular shape (e.g. a polygon) directly.
(Although note that support function methods can handle polygons in vertex representation
since the support function is easy to compute).

```@example linear_methods
A = [0 1; -1 0.]

X0 = VPolygon([[1.0, 0.0], [1.2, 0.0], [1.1, 0.2]])
prob = @ivp(x' = A*x, x(0) ∈ X0)
sol = solve(prob, tspan=(0.0, 2.0), alg=VREP(δ=1e-2))

plot(sol, vars=(1, 2), ratio=1., xlab="x", ylab="y")
plot!(X0, vars=(1, 2), color=:red)
```

!!! tip "Performance tip"
    Instead of usual Julia arrays, it is possible to use heap-allocated statically
    sized arrays with this algorithm; it suffices to set the flag `static=true`, e.g.
    `alg=VREP(δ=1e-2, static=true)`. For reference, the output of `@btime`
    for the flowpipe portion computed above runs in `189.297 μs (2890 allocations: 261.77 KiB)`
    compared to  `78.845 μs (1704 allocations: 191.30 KiB)` if we use static arrays.

In order to use `VREP` in dimensions higher than two, you have to specify a polyhedral
backend. See the documentation of [`VREP`](@ref) for details.

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

## References
