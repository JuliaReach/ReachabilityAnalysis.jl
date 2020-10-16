# Linear ordinary differential equations

In this section we show how to solve initial value problems with sets of initial
conditions. First we consider the simple *scalar* equation $x'(t) = -x(t)$,
where $x(0)$ may be any point in a given initial interval $\mathcal{X}_0$.
Then we consider the [spring-mass system](https://en.wikipedia.org/wiki/Simple_harmonic_motion#Examples),
a second order linear ODE studied in introductory physics courses. For that system
we show an example of how to compute and project the flowpipe, and then plot the
variables of interest.

## Example

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

This example illustrates that the solution process mainly consists of three steps:

(1) **Formulating** the mathematical problem, in the form of an initial-value problem
    with possibly uncertain initial states or inputs.

(2) **Solving** the problem, either with the default algorithm or specifying the algorithm
    and some of its options.

(3) **Extracting** the results, either to visualize with a plot, or to project onto
    the relevant variables for further study.

Below we give some further details of this solution step for the simple scalar equation
presented above.

### Problem formulation

Solution process

To illustrate the solution process, we consider a spring-mass system illustrated in the following figure.

The mass is $m$ and the elastic constant of the spring is $k$.

### Solving the initial-value problem

In this case the system is not given as a set of first-order ODEs, so we will make that transformation as a first step.

Transforming higher-order into a first-order system.

Formulating the mathematical problem involves writing the system as a first-order

### Analyzing the solution

## Spring-mass system

In the following example we consider a spring-mass system which is a linear ODE
with two-degrees of freedom.

## Oscillatory systems

Second order systems of the form
```math
    Mx''(t) + Cx'(t) + Kx(t) = f(t)
```
can be automatically transformed to linear systems and solved using linear reachability
solvers. For example, let's solve the damped oscillating system without a forcing
term,
```math
    x''(t) + 0.5~x'(t) + 4x(t) = 0, \qquad x(0) ∈ [0.7 .. 1.3] × [0.7 .. 1.3]
```

```@example second_order_damped
using ReachabilityAnalysis, Plots

# x'' + 0.5x' + 4x = 0
sys = SecondOrderLinearContinuousSystem(hcat([1.0]), hcat([0.5]), hcat([4.0]))

B0 = BallInf(ones(2), 0.3)
sol = solve(@ivp(sys, x(0) ∈ B0), tspan=(0.0, 10.0), alg=GLGM06(δ=0.01));

plot(sol, vars=(0, 1), lw=.2, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), lw=.2, xlab="time", lab="v(t)")
```
