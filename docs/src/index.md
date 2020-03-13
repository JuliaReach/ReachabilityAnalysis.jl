```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# ReachabilityAnalysis.jl

*Methods to compute sets of states reachable by dynamical systems.*

!!! note
    `ReachabilityAnalysis` is still under development. If you have questions,
    find a bug, or have ideas for improvements, feel free to open an issue or make
    a pull request on the [project's GitHub page](https://github.com/mforets/ReachabilityAnalysis.jl).
    You can also find us at [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).

## Introduction

This library implements reachability analysis methods for systems of ordinary
differential equations (ODEs), for both continuous and hybrid dynamical systems.

In simple terms, reachability analysis is concerned with studying the sets of states
that a system can reach, starting from a set of initial states and under the
influence of a set of input trajectories and parameter values.
More technically, we define the reachable set at a given time point
$\delta \in \mathbb{R}$, also known as the *reach-set* for the ODE
```math
x' = f(x(t), u(t), p(t), t),
```
as given by
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
\mathcal{F}(0 .. δ) = ⋃_{t \in [0, δ]} \mathcal{R}(δ).
```
Reachability methods are used to compute rigorous approximations of the flowpipe
for continuous or hybrid systems, in bounded time or unbounded time horizon.
Here we use the term *rigorous* in the formal, or mathematical sense, that no
solution "escapes" the flowpipe, for any trajectory that satisfies the constraints
(initial states, inputs, and noise).

On the other hand, the amount of computation required depends heavily on the
particular problem statement. One notable example is *safety verification*,
which simply stated requires to prove that the flowpipe does not intersect a region
of "bad states". In this setting, one can often reason about the flowpipe lazily,
i.e. without actually computing it in full.

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

## Features

The following types of ODEs are currently supported:

- Continuous ODEs with linear dynamics
- Continuous ODEs with non-linear dynamics
- Hybrid systems with piecewise-affine dynamics
- Hybrid systems with non-linear dynamics

## Application domains

Reachability analysis has applications in diverse domains such as:

- **Safety verification:** Determining whether a system is safe, i.e. if it does not
  enter into a region of unsafe sets. Typical applications are assessing the critical
  distance between autonomous vehicles or robots, or critical concentration of
  chemicals in a reactor.

- **Validation of control strategies:** Checking if the system trajectories stay in a
  region around a reference trajectory, or reach a goal region around a setpoint,
  for any admissible value of the non-deterministic inputs, initial conditions
  or noise.

- **Controller synthesis:** Finding parameter sets of controllers tha satisfy
  safety or performance constraints.


- **Deep neural network verification:** Providing formal guarantees for the network
  behavior subject to perturbations in the inputs, e.g. detecting that small
  changes in an input image do not cause the network o misclassify it.

We refer to the technical literature for further applications.
