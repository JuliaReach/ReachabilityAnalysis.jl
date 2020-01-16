# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/mforets/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/mforets/ReachabilityAnalysis.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/mforets/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/mforets/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/mforets/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/mforets/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

A Reachability Analysis suite for Dynamical Systems.

## Resources

- [Manual](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
- [Contributing](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/about/#Contributing-1)
- [Release notes of tagged versions](https://github.com/mforets/ReachabilityAnalysis.jl/releases)
- [Release log of the development version](https://github.com/JuliaReach/ReachabilityAnalysis.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/publications/)
- [Citations](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/citations/)

## Features

`ReachabilityAnalysis.jl` can handle ordinary differential equations (ODEs) with uncertain
initial states, uncertain parameters and non-deterministic inputs.

The following types of systems are supported:

- Continuous ODEs with linear dynamics
- Continuous ODEs with non-linear dynamics
- Hybrid systems with piecweise-affine dynamics
- Hybrid systems with non-linear dynamics

We suggest to explore the Examples section of the online documentation for examples illustrating how to use the library.

Full references to the scientific papers presenting the algorithms implemented in this package can be generally found in the source code as well as in the [online documentation](http://juliareach.github.io/ReachabilityAnalysis.jl/dev/).

## Installation

### Compatibility

This package is compatible with the current stable release of Julia v1.3, with
the long-term release v1.0.5 and with all intermediate versions. The package is expected to work
under Linux OS, MacOS and Windows. If you find problems installing,
do not hesitate to open an issue in the [issue tracker](https://github.com/mforets/ReachabilityAnalysis.jl/issues).

Refer to the [official documentation](https://julialang.org/downloads) on how to
install and run Julia in your system.

>  :book: This package is still a work-in-progress. It grew out of the package
    [JuliaReach/Reachability.jl](https://github.com/JuliaReach/Reachability.jl).
    If you are interested to know more, feel free to join the chat room
    [JuliaReach gitter channel](https://gitter.im/JuliaReach/Lobby).

### Installation

Depending on your needs, choose an appropriate command from the following list
and enter it in Julia's REPL. To activate the `pkg` mode, type `]` (and to leave it, type `<backspace>`).

#### [Install the latest release version](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-registered-packages-1) (recommended for first-time use)

```julia
pkg> add https://github.com/mforets/ReachabilityAnalysis.jl
```

#### Install the latest development version

```julia
pkg> add https://github.com/mforets/ReachabilityAnalysis.jl#master
```

#### [Clone the package for development](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1)

```julia
pkg> dev https://github.com/mforets/ReachabilityAnalysis.jl
```

### Dependencies

The set-based computations are handled by the library [`LazySets.jl`](https://github.com/JuliaReach/LazySets.jl),
which is part of the [JuliaReach](https://github.com/JuliaReach/) framework.

The full list of dependencies can be found in the `Project.toml` file in this repository,
in the `[deps]` section. These dependencies are installed automatically when
you `add` `ReachabilityAnaysis.jl` as explained in the previous paragraphs.
