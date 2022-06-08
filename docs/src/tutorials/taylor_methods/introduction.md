```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Taylor methods

## Univariate example

We begin the section with a small example that mixes intervals, the sine function, and differential equations. It is easy to check that the function

```math
x(t) = x_0 \exp\left(\cos(t) - 1\right),\qquad x_0 \in \mathbb{R},
```
solves the following differential equation:

```math
x'(t) = -x(t) ~ \sin(t),\qquad t \geq 0.
```

Standard integration schemes fail to produce helpful solutions if the initial state is an interval. We illustrate this point
by solving the given differential equation with the `Tsit5` algorithm from [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) suite.

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
There is no plot recipe readily available so we create it by hand using [`LazySets.jl`](https://github.com/JuliaReach/LazySets.jl).

```@example nonlinear_univariate
using LazySets, Plots
using LazySets: Interval

out = [Interval(sol.t[i]) × Interval(sol.u[i][1]) for i in 1:20]

plot(out, xlab="t", ylab="x(t)", lw=3.0, alpha=1., c=:black, marker=:none, lab="", title="Standard integrator with an interval initial condition")
```
Such divergence observed in the solution is due to using an algorithm which doesn't take into account the dependency problem when the initial condition is an interval.

On the contrary, specialized algorithms can handle this case without noticeable wrapping effect, producing a sequence of sets whose union covers the true solution for all initial points. The corresponding solution using [`ReachabilityAnalysis.jl`](https://github.com/JuliaReach/ReachabilityAnalysis.jl) can be obtained as follows:

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

Several methods exist to analyze the solution. For example, the number of reach-sets can be obtained using
`length` (or [`numrsets`](@ref) more generally for hybrid solutions):

```@example nonlinear_univariate
# number of computed reach-sets
length(sol)
```
The type of reach-set representation as well as the set representation can be obtained as well:
```@example nonlinear_univariate
# type of the reach-set representation
rsetrep(sol)
```
```@example nonlinear_univariate
# type of the set representation
setrep(sol)
```

Solutions implement the array interface hence it is possible to index them: 
```@example nonlinear_univariate
# first reach-set of the flowpipe
sol[1]
```
The associated time span can be obtained with [`tspan`](@ref):
```@example nonlinear_univariate
# time span of this reach-set
sol[1]
```
Similarly we can print the last reach-set: 
```@example nonlinear_univariate
sol[end-1:end]
```

```@example nonlinear_univariate
tspan(sol[end])
[14.9499, 15]
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

Finally, note that the solution structure, apart from storing the resulting flowpipe, contains information
about the algorithm that was used to obtain it:

```@example nonlinear_univariate
sol.alg
```
