# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/mforets/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/mforets/ReachabilityAnalysis.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/mforets/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/mforets/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/mforets/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

>  :book: This package is still a work-in-progress. It grew out of the package
    [JuliaReach/Reachability.jl](https://github.com/JuliaReach/Reachability.jl).
    If you are interested to know more about the project, feel free to open an issue or
    join the chat room [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby).

NOTE: i'm currently using `mforets/dev` branch in LazySets.jl.

## Resources

- [Manual](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
- [Contributing](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/about/#Contributing-1)
- [Release notes of tagged versions](https://github.com/mforets/ReachabilityAnalysis.jl/releases)
- [Release log of the development version](https://github.com/JuliaReach/ReachabilityAnalysis.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/publications/)
- [Citations](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/citations/)

## Features

Reachability Analysis is concerned with computing rigorous approximations of the set
of states reachable by a dynamical system. The scope of the package are systems
modeled by ordinary differential equations (ODEs) with uncertain initial states,
uncertain parameters or non-deterministic inputs.

The following types of systems are supported:

- Continuous ODEs with linear dynamics
- Continuous ODEs with non-linear dynamics
- Hybrid systems with piecweise-affine dynamics
- Hybrid systems with non-linear dynamics

We suggest to explore the Examples section of the online documentation for examples
illustrating how to use the library.

References to the scientific papers presenting the algorithms implemented in this
package can be found in the source code and in
the [online documentation](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/).

## Installation

Once you have installed Julia in your system, open a Julia session, activate the
`pkg` mode (remember that to activate the `pkg` mode in Julia's REPL, you need to type `]`,
and to leave it, type `<backspace>`), and enter:

```julia
pkg> add ReachabilityAnalysis.jl
```

## Examples

```julia
using ReachabilityAnalysis

prob = @ivp(x' = 1.01x, x(0) ∈ 0 .. 1)
sol = solve(prob, T=1.0, GLGM06())

using Plots

plot(sol)
```

```julia
using ReachabilityAnalysis

A = [1 0; 0 -1]
B = [1, 1]
U = Interval(-0.1, 0.1)
prob = @ivp(x' = Ax + Bu, u ∈ U, x(0) ∈ X0)
sol = solve(prob, T=1.0, GLGM06())

using Plots

plot(sol)
```
