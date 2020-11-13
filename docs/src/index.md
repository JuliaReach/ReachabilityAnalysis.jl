```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# ReachabilityAnalysis.jl

*Methods to compute sets of states reachable by dynamical systems.*

!!! note
    `ReachabilityAnalysis` is still under development. If you have questions,
    find a bug, or have ideas for improvements, feel free to open an issue or make
    a pull request on the [project's GitHub page](https://github.com/mforets/ReachabilityAnalysis.jl).
    You can also find us at the [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).

## Features

This library implements reachability analysis methods for systems of ordinary
differential equations (ODEs), for both continuous and hybrid dynamical systems.

The following types of ODEs are currently supported:

- Continuous ODEs with linear dynamics.
- Continuous ODEs with non-linear dynamics.
- Hybrid systems with piecewise-affine dynamics.
- Hybrid systems with non-linear dynamics.

For hybrid systems, the transitions may be space-triggered or time-triggered
(or both).

In all the problems mentioned above, the library can handle uncertainties in
the sets of initial states, inputs, or parameter variation.

We refer to the `Algorithms` section for detailed descriptions of the algorithms
available, as well as the references to the technical literature.

## What can `ReachabilityAnalysis` do?

- Solve the same ODE repeatedly for an *infinite* number of different initial conditions.

## How to install the package?

You can easily install `ReachabilityAnalysis.jl` from your Julia console
(press `]` in the Julia REPL to enter the `pkg>` mode):

```julia
pkg> add ReachabilityAnalysis
```

## Application domains

Reachability analysis has applications in diverse domains such as:

- **Safety verification:** Determining whether a system is safe, i.e. to verify
  if it does not enter into a region of unsafe sets. Typical applications are
  assessing the critical distance between autonomous vehicles or robots, or critical
  concentration of chemicals in a reactor.

- **Validation of control strategies:** Checking if the system trajectories stay in a
  region around a reference trajectory, or reach a goal region around a setpoint,
  for any admissible value of the non-deterministic inputs, initial conditions
  or noise.

- **Controller synthesis:** Finding parameter sets of controllers that satisfy
  safety or performance constraints.

- **Deep neural network verification:** Providing formal guarantees for the network
  behavior subject to perturbations in the inputs, e.g. detecting that small
  changes in an input image do not cause the network o misclassify it.

We refer to the technical literature for further applications.

## What is the verification problem?


## A warning note

Users of the library should have in mind that parameter tuning is an
essential ingredient to the successful application of reachability analysis.
This is in part due to the "cutting-edge research" aspect
of the methods available, i.e. the methods haven't yet been "battle-tested".
Moreover, finding good algorithm heuristics -- specially for hybrid systems -- is a hard
problem by itself.

On the other hand, our goal in designing and building the tools around `JuliaReach` has
been to make the default settings already work, or at least, work reasonably well across
a broad range of models: linear, hybrid, nonlinear, and so on. In the future, we hope that
users won't have to know the details of how the different set-based reachability methods
work under-the-hood in order to get new results. Eventually, users should know how
to adjust *some* parameters associated to the algorithms in order to improve the precision,
to formulate and solve verification problems, etc. inasmuch as effectively applying
standard numerical integration requires some basic knowledge about stepsize control,
interpolation, error calculation, etc.


## How to read this manual

Please note that the manual contains a [Frequently Asked Questions (FAQ)](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/man/faq/) section with a selection of different questions and answers, many of them
with working code so that you can also experiment with the answers.
