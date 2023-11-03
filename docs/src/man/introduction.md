## Introduction

In simple terms, reachability analysis is concerned with studying the sets of states
that a system can reach, starting from a set of initial states and under the
influence of a set of input trajectories and parameter values.

Our first example is the one-dimensional scalar differential equation

```math
x'(t) = -x(t),\qquad 0 ≤ t ≤ T = 4.0,
```
with initial condition $x(0) ∈ [0.45, 0.55]$. We can solve the problem as follows:

```@example linear_scalar
using ReachabilityAnalysis, Plots

# define the initial-value problem
prob = @ivp(x' = -x, x(0) ∈ 0.45 .. 0.55)

# solve it
sol = solve(prob, T=4.0)

# plot the solution, where the index 0 corresponds to the "time" variable
plot(sol, vars=(0, 1), label="Flowpipe", xlab="t", ylab="x(t)", linewidth=0.3)

# also plot the analytic solution
f(t, x0) = x0 * exp(-t)
plot!(t -> f(t, 0.45), xlims=(0, 4), label="Analytic sol., x(0) = 0.45", color="red")
plot!(t -> f(t, 0.55), xlims=(0, 4), label="Analytic sol., x(0) = 0.55", color="red")
```
In practice, analytic solutions of ODEs are unknown. However, in this simple case
we knew that for an initial point $x_0 \in \mathbb{R}$, the solution is
$x(t) = x_0 e^{-t}$ so we plotted the trajectories associated to the extremal
values in the given initial interval.

We have seen how to compute the set of all reachable stataes within
the time-span $0 ≤ t ≤ T = 4.0$. The reachability approach guarantees that
the trajectories are always included in the resulting *flowpipe*. In the following
section we introduce the terminology used in reachability analysis and formulate the
reachability problem for a general class of mathematical systems. Further examples on
how to use the library for other types of systems can be found in later sections of the manual.

## Terminology

The reachable set at a given time point $\delta \in \mathbb{R}$, also known as the
*reach-set* for the ODE
```math
x'(t) = f(x(t), u(t), p(t), t),
```
is defined by
```math
\mathcal{R}(δ) := \left\{ x(δ) = \int_0^δ f(x(t), u(t), p(t), t) dt, x(0) ∈ X_0, u(t) ∈ \mathcal{U}, p(t) ∈ \mathcal{P} \right\}.
```
Here $X_0$ denotes the set of initial states, $\mathcal{U}$ denotes the input set,
and $\mathcal{P}$ denotes the parameter values. For practical problems, the set
$\mathcal{R}(δ)$ cannot be obtained exactly, and reachability methods aim at
computing suitable over-approximations (or under-approximations) of it.

We define the reachable set associated to a time interval $[0, δ]$,
also known as the *flowpipe*, as
```math
\mathcal{F}([0, δ]) = ⋃_{t \in [0, δ]} \mathcal{R}(t).
```
Reachability methods are used to compute rigorous approximations of the flowpipe
for continuous or hybrid systems, in bounded time or unbounded time horizon.
Here we use the term *rigorous* in the formal, or mathematical sense, that no
solution "escapes" the flowpipe, for any trajectory that satisfies the constraints
(initial states, inputs, and noise).

## Safety verification

On the other hand, the amount of computation required depends heavily on the
particular problem statement. One notable example is *safety verification*,
which simply stated requires to prove that the flowpipe does not intersect a region
of "bad states". In this setting, one can often reason about the flowpipe lazily,
i.e. without actually computing it in full.

## Hybrid systems

Up to now we have discussed about the continuous case only, but there is a rich
literature in hybrid systems reachability; *hybrid* here means those dynamical
systems which are given by one or more continuous-time dynamics (often, systems
of ODEs in each mode or location) coupled with discrete transitions between
continuous modes. In our context it is standard to model these systems using the
terminology of *hybrid automata*, and we also model hybrid systems with such framework
in this library. The concept of reach-set, flowpipe and safety verification are
naturally extended to hybrid automata, although there is the additional complication
that the flowpipe must include the behaviors for all possible transitions between
discrete modes that are compatible with the dynamics.
