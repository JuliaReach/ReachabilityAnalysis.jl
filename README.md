# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/mforets/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/mforets/ReachabilityAnalysis.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://mforets.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/mforets/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/mforets/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/mforets/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


Reachability Analysis is concerned with computing rigorous approximations of the set
of states reachable by a dynamical system. The scope of the package are systems
modeled by ordinary differential equations (ODEs) with uncertain initial states,
uncertain parameters or non-deterministic inputs.

## Resources

- [Reference manual](https://mforets.github.io/ReachabilityAnalysis.jl/dev/)
- [Citation](#citation)

## Features

The following types of systems are supported.

### Implemented

- Continuous ODEs with linear dynamics (GLGM06, INT) :heavy_check_mark:
- Continuous ODEs with linear dynamics and parametric uncertainty (ASB07) :heavy_check_mark:
- Continuous ODEs with non-linear dynamics (TMJets) :heavy_check_mark:
- Hybrid systems with clocked linear dynamics (HACLD1 with GLGM06 or ASB07) :heavy_check_mark:

### Planned

- Continuous ODEs with linear dynamics (BFFPSV18)
- Hybrid systems with piecweise-affine dynamics
- Hybrid systems with non-linear dynamics

Research credit and full references to the scientific papers presenting the algorithms
implemented in this package can be found in the source code for each algorithm and in the
[References](https://mforets.github.io/ReachabilityAnalysis.jl/dev/references/) section of the online documentation.

## Installation

Open a Julia session and activate the
`pkg` mode (to activate the `pkg` mode in Julia's REPL, type `]`,
and to leave it, type `<backspace>`), and enter:

```julia
pkg> add ReachabilityAnalysis
```

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
