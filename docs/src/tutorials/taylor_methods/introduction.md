```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Introduction

This section provides a quickstart to the Taylor-models reachability method with a *univariate oscillator* as running example.
It is a small example that mixes intervals, the sine function, and differential equations. It is easy to check that the function

```math
x(t) = x_0 \exp\left(\cos(t) - 1\right),\qquad x_0 ∈ \mathbb{R},
```
solves the following differential equation:

```math
x'(t) = -x(t) ~ \sin(t),\qquad t ≥ 0.
```

Standard integration schemes fail to produce helpful solutions if the initial state is an interval. We illustrate this point
by solving the given differential equation with the `Tsit5` algorithm from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) suite.

```@example nonlinear_univariate
using OrdinaryDiffEq, IntervalArithmetic

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
There is no plot recipe readily available so we create it by hand using [`LazySets.jl`](https://github.com/JuliaReach/LazySets.jl).

```@example nonlinear_univariate
using LazySets, Plots
using LazySets: Interval

out = [Interval(sol.t[i]) × Interval(sol.u[i][1]) for i in 1:20]

fig = plot(out, xlab="t", ylab="x(t)", lw=3.0, alpha=1., c=:black, marker=:none, lab="", title="Standard integrator with an interval initial condition")

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```
The divergence observed in the solution is due to using an algorithm which doesn't specialize for intervals hence suffers from dependency problems.

However, specialized algorithms can handle this case without noticeable wrapping effect, producing a sequence of sets whose union covers the true solution for all initial points. We use the [`ReachabilityAnalysis.jl`](https://github.com/JuliaReach/ReachabilityAnalysis.jl) interface to the algorithm [`TMJets21a`](@ref), which is itself implemented in
[`TaylorModels.jl`](https://github.com/JuliaIntervals/TaylorModels.jl).

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
fig = plot(sol, vars=(0, 1), xlab="t", ylab="x(t)", title="Specialized (Taylor-model based) integrator")

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

It is illustrative to plot the computed flowpipe and the known analytic solution
for a range of initial conditions that cover the range of the initial interval `X0 = -1 .. 1`.
We can do so by overlaying exact solutions taken uniformly from the initial interval.

```@example nonlinear_univariate
analytic_sol(x0) = t -> x0 * exp(cos(t) - 1.0)

dt = range(0, 15, length=100)
x0vals = range(-1, 1, length=25)

fig = plot()
plot!(fig, sol, lw=0.0, vars=(0, 1), xlab="t", ylab="x(t)", title="Specialized (Taylor-model based) integrator")
[plot!(fig, dt, analytic_sol(x0).(dt), lab="", c=:magenta, lw=2.0, xlab="t", ylab="x(t)") for x0 in x0vals]
fig

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

[`ReachabilityAnalysis.jl`](https://github.com/JuliaReach/ReachabilityAnalysis.jl) defines several methods to analyze the solution.
For example, the number of reach-sets can be obtained using `length` (more generally, use `numrsets` for hybrid systems):

```@example nonlinear_univariate
# number of computed reach-sets
length(sol)
```

Solutions implement the array interface hence it is possible to index them:
```@example nonlinear_univariate
# first reach-set of the flowpipe
sol[1]
```
Note that `sol[1]` is an expansion in *time* whose coefficients are polynomials in *space*.
The associated time span can be obtained with [`tspan`](@ref):
```@example nonlinear_univariate
# time span of this reach-set
tspan(sol[1])
```
Similarly we can print the last reach-set:
```@example nonlinear_univariate
sol[end]
```

```@example nonlinear_univariate
tspan(sol[end])
```
It is also possible to filter solutions by the time variable using parentheses:

```@example nonlinear_univariate
# find the reach-set(s) that contain 1.0 in their time span
aux = sol(1.0)
```

```@example nonlinear_univariate
tspan(aux)
```

Filtering by the time span also works with time intervals; in the following example, we obtain all reach-sets whose time span with the interval
`1.0 .. 3.0` is non-empty:

```@example nonlinear_univariate
length(sol(1.0 .. 3.0))
```

On the other hand, evaluating over a given time point or time interval can be achieved using the `evaluate` function:

```@example nonlinear_univariate
# evaluate over the full time interval ≈ [0, 0.100777]
evaluate(sol[1], tspan(sol[1]))
```
```@example nonlinear_univariate
# evaluate over the final time point 0.10077670723811318
evaluate(sol[1], tend(sol[1]))
```
```@example nonlinear_univariate
# evaluate the final reach-set at the final time point
evaluate(sol[end], 15.0)
```

Finally, note that the solution structure, apart from storing the resulting flowpipe, contains information
about the algorithm that was used to obtain it:

```@example nonlinear_univariate
sol.alg
```

The type of reach-set representation as well as the set representation can be obtained using `rsetrep` and [`setrep`](@ref) respectively.
```@example nonlinear_univariate
# type of the reach-set representation
rsetrep(sol)
```
```@example nonlinear_univariate
# type of the set representation
setrep(sol)
```
The following section discusses more details about those representations.
