# Linear ordinary differential equations

## Introduction

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
x0sample = [0.45, 0.50, 0.55]
[plot!(fig, t -> x0*exp(-t), xlims=(0, 4), label="Analytic with x(0) = $x0") for x0 in x0sample]
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

---

This section is about initial-value problems for linear differential equations of the form

```math
x'(t) = Ax(t) + u(t), \qquad x(0) \in \mathcal{X}_0, u(t) \in \mathcal{U}(t),
```
where ``x(t) \in \mathbb{R}^n`` is the state vector, ``A \in \mathbb{R}^{n \times n}``
is the state (or coefficients) matrix, ``\mathcal{X}_0`` is the set of initial states,
and ``\mathcal{U}(t)`` is the (possibly time-varying) set of uncertain inputs.
Although many systems can be formulated in this manner, arguably most systems of interest
in science and engineering have additional complexities, e.g. the state matrix `A`
can depend on time, the inputs can be state-dependent, or more generally right-hand
side can contain powers or other nonlinear functions of the state `x`, etc. Moreover,
in several situations of interest the coefficients of ``A`` are only known approximately.
On the other hand, studying linear systems is of central interest for our purposes

!!! note "Extension tip"
    Please note that equations of the form ``x'(t) = Ax(t) + Bu(t)`` are also considered
    in this context, since . . . .

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


## Zonotope methods

## Support function methods

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

## Template polyhedra 



Further details on the methods presented in this section is presented in the section
[Exploiting structure](@ref).

## Working with solutions

Once the flowpipe has been computed, we are interested in extracting information
from it.


In the following example we consider a spring-mass system which is a linear ODE
with two-degrees of freedom.

Analyzing the solution

## Oscillating systems

Second order linear systems of the form
```math
    Mx''(t) + Cx'(t) + Kx(t) = f(t),
```
assuming $M$ is invertible, can be solved using linear reachability
methods.

### Free oscillations

For example, let's solve the damped oscillating system without a forcing
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
Where we have chosen the zonotope-based algorithm `GLGM06` with step-size`δ=0.01` and plotted the flowpipe
for $x(t)$ and $v(t)$ variables.

### Constant forcing term

For example, let's solve the damped oscillating system without a forcing
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
Where we have chosen the zonotope-based algorithm `GLGM06` with step-size`δ=0.01` and plotted the flowpipe
for $x(t)$ and $v(t)$ variables.

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

```@example
F = 3.0
A = [ 0    1.0     0   0;
     -4    -3.5    1   0;
      0     0      0   F;
      0     0     -1   0]
sys = @system(x' = A * x)

B0 = BallInf(ones(2), 0.3) × Singleton([0.0]) × Interval(0.8, 2.1)
sol = solve(@ivp(sys, x(0) ∈ B0), tspan=(0.0, 15.0), alg=GLGM06(δ=0.01));

plot(sol, vars=(0, 1), c=:magenta, lw=.0, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), c=:green, lw=.0, xlab="time", lab="v(t)", legend=:bottomright)
```
