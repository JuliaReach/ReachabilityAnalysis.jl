# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/JuliaReach/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/ReachabilityAnalysis.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/juliareach/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/juliareach/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/juliareach/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## ✨  What is Reachability Analysis?

Reachability analysis is concerned with computing rigorous approximations of the set
of states reachable by a dynamical system. In the scope of this package are systems
modeled by ordinary differential equations (ODEs) with uncertain initial states,
uncertain parameters or non-deterministic inputs. The package also considers the
extension to so-called hybrid systems where the dynamics changes with discrete events.

## 🎯  Table of Contents

* [Installation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--installation)
* [Documentation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--documentation)
* [Features](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--features)
* [Examples Gallery](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--examples)
* [References](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--references)

## 💾  Installation

Open a Julia session and activate the
`pkg` mode (to activate the `pkg` mode in Julia's REPL, type `]`,
and to leave it, type `<backspace>`), and enter:

```julia
pkg> add ReachabilityAnalysis
```

## 📙  Documentation

See the [Reference Manual](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/) for introductory material, examples and API reference.

📌 **Need help?** Have any question, or wish to suggest or ask for a missing feature?
Do not hesitate to open a [**issue**](https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues) or join the JuliaReach gitter chat: [![join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge), or send an [email](mailto:mforets@gmail.com) to the developers.


## 🎨  Features

The following types of systems are supported.

- Continuous ODEs with linear dynamics :heavy_check_mark:
- Continuous ODEs with linear dynamics and parametric uncertainty :heavy_check_mark:
- Continuous ODEs with non-linear dynamics  :heavy_check_mark:
- Hybrid systems with piecewise-affine dynamics :heavy_check_mark:
- Hybrid systems with non-linear dynamics :heavy_check_mark:
- Hybrid systems with clocked linear dynamics :heavy_check_mark:


## 🐾  Examples

|         |      |
|:--------:|:-----:|
| ![quadrotor](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-Quadrotor.png?raw=true) [Quadrotor altitude control](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/quadrotor/) |  ![LVHybrid](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-LotkaVolterra.png?raw=true) Lotka-Volterra with tangential guard crossing|
|         |      |
| ![LaubLoomis](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-LaubLoomis.png?raw=true) [Laub-Loomis model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/laub_loomis/)    | ![PD](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-ProductionDestruction.png?raw=true)<br> Production-Destruction model|
|         |      |
|![vanderpol](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-VanDerPol.png?raw=true) Coupled van der pol oscillator model  | ![spaccecraft](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-Spacecraft.png?raw=true) [Spacecraft rendez-vous model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/spacecraft/) |


## 📜  References

Research credit and full references to the scientific papers presenting the algorithms
implemented in this package can be found in the source code for each algorithm and in the
[References](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/references/) section of the online documentation.

The following articles have used `ReachabilityAnalysis.jl` for their research:

- *Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.* Marcelo Forets, Daniel Freire, Christian Schilling, 2020. [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325).

---

If you use this package for your research, we kindly ask you to cite the following paper, see [CITATION.bib](http://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/CITATION.bib). Moreover, please also cite the appropriate originl references to the algorithm(s) used.



