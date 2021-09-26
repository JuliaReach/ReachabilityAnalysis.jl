```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Taylor models methods

## Univariate example

We begin the section with a small example that mixes intervals, the sine function, and differential equations. It is easy to check that the function

```math
x(t) = x_0 \exp\left(\cos(t) - 1\right),\qquad x_0 \in \mathbb{R},
```
solves the following differential equation:

```math
x'(t) = -x(t) ~ \sin(t),\qquad t \geq 0.
```

Standard integration schemes fail to produce helpful solutions if the initial state is an interval,
for different reason, e.g. dependency problems.

```@example nonlinear_univariate
using DifferentialEquations, IntervalArithmetic

# initial condition
x₀ = [-1 .. 1]

# define the problem
function f(dx, x, p, t)
    dx[1] = -x[1] * sin(t)
end

# pass to solvers
prob = ODEProblem(f, x₀, (0.0, 2.0))
sol = solve(prob, Tsit5(), adaptive=false, dt=0.05, reltol=1e-6)
nothing # hide
```
There is no plot recipe readily available so we create it by hand using LazySets.

```@example nonlinear_univariate
using LazySets, Plots
using LazySets: Interval

out = [Interval(sol.t[i]) × Interval(sol.u[i][1]) for i in 1:20]

plot(out, xlab="t", ylab="x(t)", lw=3.0, alpha=1., c=:black, marker=:none, lab="", title="Standard integrator with an interval initial condition")
```

On the contrary, specialized algorithms can handle this case without noticeable wrapping effect, producing a sequence of sets whose union covers the true solution for all initial points.
The corresponding solution using reachability analysis is simple; .

```@example nonlinear_univariate
using ReachabilityAnalysis, Plots

# define the model (same as before)
function f(dx, x, p, t)
    dx[1] = -x[1] * sin(1.0 * t)
end

# define the set of initial states
X0 = -1 .. 1

# define the initial-value problem
prob = @ivp(x' = f(x), x(0) ∈ X0, dim=1)

# solve it using a Taylor model method
sol = solve(prob, alg=TMJets21a(abstol=1e-10), T=15)

# visualize the solution in time
plot(sol, vars=(0, 1), xlab="t", ylab="x(t)", title="Specialized (Taylor-model based) integrator")
```

It is illustrative to plot the computed flowpipe and the known analytic solution
for a range of initial conditions that cover the range of `U0`. We can do so by
overlaying exact solutions taken uniformly from the initial interval.

```@example nonlinear_univariate
analytic_sol(x0) = t -> x0 * exp(cos(t) - 1.0)

dt = range(0, 15, length=100)
x0vals = range(-1, 1, length=25)

fig = plot()
plot!(fig, sol, lw=0.0, vars=(0, 1), xlab="t", ylab="x(t)", title="Specialized (Taylor-model based) integrator")
[plot!(fig, dt, analytic_sol(x0).(dt), lab="", c=:magenta, lw=2.0, xlab="t", ylab="x(t)") for x0 in x0vals]
fig
```

```@example nonlinear_univariate
# number of computed reach-sets
length(sol)
```

```@example nonlinear_univariate
# reach-set representation
typeof(sol[1])
```

```@example nonlinear_univariate
set(sol[1])
```

```@example nonlinear_univariate
sol_univariate = deepcopy(sol); # for further use in this notebook
```


## Overview

Reachability analysis applies to nonlinear systems, i.e. those where the right-hand
side of the ODE is a nonlinear function of the state variables. Such systems
play a central role in applied mathematics. In this section we explain how to solve
nonlinear reachability problems using `ReachabilityAnalysis.jl` and comment on some
noteworthy differences between the user interface of nonlinear vs. linear systems.
Since nonlinear reachability methods suffer from wrapping effects, we explain some
common techniques to improve error bounds, such as splitting and refinement.
Computing accurately and efficiently the sets of states reachable by nonlinear ODEs
is a hard mathematical and computational problem. The methods available in this
library are just a small portion of the active and rapidly evolving research
literature.

In this section we take, as our running example, the well-known
[Lotka-Volterra equations](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations).
The two-dimensional Lotka-Volterra system depicts the populations change of a class of predators and a class of preys. The growth rate of preys' population $x$ over time is given by

```math
\dot{x} = x\cdot (\alpha - \beta \cdot y)
```
wherein  $\alpha, \beta$ are constant parameters and $y$ is the population of predators.
This equation states that the number of preys grows exponentially without predation.
On the other hand, the population growth of predators is governed by the differential
equation

```math
\dot{y} = -y\cdot (\gamma - \delta\cdot x)
```
wherein  $\gamma, \delta$ are constant parameters.

A typical choice of parameter values is $\alpha = 1.5$, $\beta = 1$,
$\gamma = 3$ and $\delta = 1$. In the next section we consider these values, and
compute the set of states reachable varying the initial condition. After, we assume
that the paramters are only known within given intervals, and compute the flowpipe
for all possible values of the parameters and initial conditions.

## Problem formulation

nonlinear systems can be computed by stating and solving an initial-value
problem for the given , similarly to the case of linear systems, but using different algorithms

The first set
We introduce the vector $\

```@example lotka_volterra
using ReachabilityAnalysis

@taylorize function lotka_volterra!(du, u, p, t)
    local α, β, γ, δ = 1.5, 1.0, 3.0, 1.0
    x, y = u
    du[1] = x * (α - β*y)
    du[2] = -y * (γ - δ*x)
end
```

## Computing with Taylor models

The reachability settings are taken from [this resource](https://ths.rwth-aachen.de/research/projects/hypro/lotka-volterra/).

We consider the initial set  $x\in [4.8,5.2], y \in [1.8,2.2]$.

```@example lotka_volterra
X₀ = Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2])

prob = @ivp(x' = lotka_volterra!(x), dim: 2, x(0) ∈ X₀)
```

## Results

We compute the flowpipe using the TMJets algorithm for the time horizon $[0,5]$:

```@example lotka_volterra
sol = solve(prob, T=5.0)

setrep(sol)
```

We can change to the zonotopic overapproximation of the flowpipe using
the `overapproximate` function:

```@example lotka_volterra
sol = overapproximate(sol, Zonotope)

setrep(sol)
```

Finally we plot the solution in phase-space:

```@example lotka_volterra
using Plots

plot(sol, vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:lightblue, lab="Flowpipe")
plot!(X₀, color=:orange, lab="Xo")
```
