```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Hybrid systems

## Introduction

!!! note TO-DO
    Cleanup and complete the section.

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

Our running example is the *bouncing ball* model; although it is a very hybrid automaton,
it can be used to introduce the main notions involved in hybrid systems reachability.

![Hybrid automaton of the bouncing ball model](../assets/bouncing_ball_annotations.png)

## Formalism

TODO: definition of HA

## Hybrid solver algorithm

TODO: high-level description of the hybrid solve

## Customizing the solver

TODO: further examples on how to customize the discrete post-operator
