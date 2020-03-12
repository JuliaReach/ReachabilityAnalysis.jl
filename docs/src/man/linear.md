# Linear ordinary differential equations

In this section we begin by solving a simple scalar equation. Then we introduce
the functionality to solve systems of equations. Finally, we consider a
higher-dimensional example and some

## Example

Our first example is the one-dimensional differential equation $x'(t) = -x$ over
the time interval $0 ≤ t ≤ T$. Suppose that the initial state can be any point in
the interval $x(0) ∈ [0, 1]$.

```julia
using ReachabilityAnalysis
P = @ivp x' = -x, x(0) ∈ 0..1
sol = solve(P)

using Plots
plot(sol, label="Flowpipe")
```

The analytic solution for an initial point $x_0$ is $x(t) = e^{-t}x_0$, which
is also plotted

```julia
plot!(, t -> -exp(-t), label="Analytic solution")
```

As we have seen in this example, the solution process mainly consists of three steps:

- formulating the mathematical problem
- solving the problem
- extracting the results, eg. plotting, or projecting

In the following example we consider a higher-dimensional ODE to illustrate in some detail these steps.

## Solution process

To illustrate the solution process, we consider a spring-mass system illustrated in the following figure.

The mass is $m$ and the elastic constant of the spring is $k$.

### Problem formulation

In this case the system is not given as a set of first-order ODEs, so we will make that transformation as a first step.

Transforming higher-order into a first-order system.

Formulating the mathematical problem involves writing the system as a first-order

### Solving the problem

### Extracting the results, eg. plotting, or projecting

### Using static arrays

## Higher-dimensional ODEs

- combination of state variables to obtain an output

### Using lazy exponentiation
