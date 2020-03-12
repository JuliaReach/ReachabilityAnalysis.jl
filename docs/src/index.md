```@meta
DocTestSetup = :(using ReachabilityAnalysis)
```

# ReachabilityAnalysis.jl

`ReachabilityAnalysis.jl` is a library to compute rigorous approximations of the
sets of states reachable by dynamical systems.  The library can handle ordinary
differential equations (ODEs) with uncertain initial states, parameters or inputs.

To illustrate, here formulate the forward-time reachability problem for a
In the technical literature, this problem is called *reachability analysis* or
*validated integration*.

```math
X = \left\{ x : \int_0^t f(x, p, t, u) dt \right\}
```

In the literature, the set $X$ is called a **flowpipe**.

## Features

The following types of ODEs are currently supported:

- Continuous ODEs with linear dynamics
- Continuous ODEs with non-linear dynamics
- Hybrid systems with piecweise-affine dynamics
- Hybrid systems with non-linear dynamics

!!! warning

    This package is still a work-in-progress and it has grown out of
    [JuliaReach/Reachability.jl](https://github.com/JuliaReach/Reachability.jl).
    If you are interested to know more about the project, if you found a bug,
    or if you want to propose an idea for improvement, feel free to open an issue or
    join the chat room [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby).

## Application domains

`ReachabilityAnalysis` is a [Julia](http://julialang.org) package for approximating the
reachable states and checking safety properties of affine systems.

Reachability analysis has applications in diverse domains such as:

- verification of deep neural networks
- algorithmic verification
- xxx
- yyy
- zzz

## Benchmarks

## Roadmap
