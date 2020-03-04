# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/mforets/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/mforets/ReachabilityAnalysis.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/mforets/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/mforets/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/mforets/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


Reachability Analysis is concerned with computing rigorous approximations of the set
of states reachable by a dynamical system. The scope of the package are systems
modeled by ordinary differential equations (ODEs) with uncertain initial states,
uncertain parameters or non-deterministic inputs.

## Resources

- [Reference manual](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
- [Citation]()

## Features

The following types of systems are supported:

- Continuous ODEs with linear dynamics
- Continuous ODEs with non-linear dynamics
- Hybrid systems with piecweise-affine dynamics
- Hybrid systems with non-linear dynamics

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


>  :book: This package is still a work-in-progress. It grew out of the package
    [JuliaReach/Reachability.jl](https://github.com/JuliaReach/Reachability.jl).
    If you are interested to know more about the project, feel free to open an issue or
    join the chat room [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby).

!! NOTE: currently using `mforets/dev` branch in LazySets.jl.

## Examples

We suggest to explore the Examples section of the online documentation for examples
illustrating how to use the library.

### Linear

```julia
using ReachabilityAnalysis

prob = @ivp(x' = 1.01x, x(0) ∈ 0 .. 1)
sol = solve(prob, T=1.0, GLGM06())

using Plots

plot(sol)
```

### Higher-dimensional linear

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

### Nonlinear


### Hybrid

## Citation

If you use this package for your research, we kindly ask you to consider citing the following paper.

```
@inproceedings{bogomolov2019juliareach,
  title={JuliaReach: a toolbox for set-based reachability},
  author={Bogomolov, Sergiy and Forets, Marcelo and Frehse, Goran and Potomkin, Kostiantyn and Schilling, Christian},
  booktitle={Proceedings of the 22nd ACM International Conference on Hybrid Systems: Computation and Control},
  pages={39--44},
  year={2019}
}
```
