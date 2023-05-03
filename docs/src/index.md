```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# ReachabilityAnalysis.jl

## What is reachability analysis?

Reachability analysis is a numerical method that aims to compute sets of states
reachable by dynamical systems from all initial states and for all admissible inputs
and parameters.

## Where is reachability being applied?

Reachability analysis is of interest in diverse domains including automotive,
aerospace, power systems, system biology, analog/mixed-signal circuits, robotics,
and safe machine learning.

Below we briefly comment on a few remarkable problems where reachability is being
applied. We refer to the technical literature for further applications;
see [Bibliography](@all_ref) for a selection of papers to get you started!
An up-to-date review with relevant references can be found in [[AFG20]](@ref).

**Formal verification.** Determining whether a system is safe, i.e. to verify
if it does not enter into a region of unsafe sets. Typical applications are
assessing the critical distance between autonomous vehicles or robots, or critical
concentration of chemicals in a reactor.

**Validation of control strategies.** Checking if the system trajectories stay in a region around a reference trajectory, or reach a goal region around a setpoint,
for any admissible value of the non-deterministic inputs, initial conditions
or noise.

**Controller synthesis.** The controller synthesis problem consists in finding parameter sets of controllers that satisfy safety or performance constraints.

**Deep neural network verification.** Providing formal guarantees for the network behavior subject to perturbations in the inputs, e.g. detecting that small changes in an input image do not cause the network o misclassify it.

## How to read this manual

Please note that the manual contains a [Frequently Asked Questions (FAQ)](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/man/faq/) section with a selection of different questions and answers, many of them
with working code so that you can also experiment with the answers.

!!! note
    If you have questions, want to chat about how to apply this package, or if you found a bug or have trouble running an example, feel free to open an issue or make a pull request on the [project's GitHub page](https://github.com/mforets/ReachabilityAnalysis.jl). You can also find us at the [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge), and at the `#reachability-analysis` stream on [JuliaLang's zulip](https://julialang.zulipchat.com) channel.

## Features

This library implements reachability analysis methods for systems of ordinary
differential equations (ODEs), for both continuous and hybrid dynamical systems.

The following types of systems are currently supported:

- Continuous ODEs with linear dynamics.
- Continuous ODEs with non-linear dynamics.
- Hybrid systems with piecewise-affine dynamics.
- Hybrid systems with non-linear dynamics.
- Hybrid systems with space-triggered or time-triggered transitions.

The library can handle uncertainties in the sets of initial states, inputs, or parameter variation.

## How to install the package?

You can easily install `ReachabilityAnalysis.jl` from your Julia console
(press `]` in the Julia REPL to enter the `pkg>` mode):

```
pkg> add ReachabilityAnalysis
```
