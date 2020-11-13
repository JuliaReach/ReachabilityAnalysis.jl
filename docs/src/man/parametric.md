```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Parametric reachability

## Exponential decay

Our first example is an initial-value problem for the one-dimensional differential
equation

```math
x'(t) = -x(t),\qquad 0 ≤ t ≤ T = 4.0,
```
with initial condition $x(0) ∈ [0.45, 0.55]$.

We can compute the flowpipe using `solve`:

```@example linear_scalar
using ReachabilityAnalysis, Plots

# define the initial-value problem
prob = @ivp(x' = -x, x(0) ∈ 0.45 .. 0.55)

# solve it
sol = solve(prob, T=4.0)

# plot the solution, where the index 0 corresponds to the "time" variable
plot(sol, vars=(0, 1), label="Flowpipe", xlab="t", ylab="x(t)", linewidth=0.3)
```

In practice, analytic solutons of ODEs are unknown. However, in this simple case
we know that for an initial point $x_0 \in \mathbb{R}$, the solution is
$x(t) = x_0 e^{-t}$. We can plot some trajectories in the same plot as the flowpipe,
to see that the trajectories are indeed inside the flowpipe, as expected.

```@example linear_scalar
f(t, x0) = x0 * exp(-t)

plot!(t -> f(t, 0.45), xlims=(0, 4), label="Analytic sol., x(0) = 0.45", color="red")
plot!(t -> f(t, 0.55), xlims=(0, 4), label="Analytic sol., x(0) = 0.55", color="red")
```

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

## Theory

Parametric reachability consist of . . .

## Example

We reconsider the example from Section (REF?), but we add an uncertain parameter $α$ that accounts for the variation in the ...


## Spring-mass system

Here we consider variations in the constant $k$ and perform reachability . . . .
