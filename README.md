# ReachabilityAnalysis.jl

[![Build Status](https://travis-ci.org/JuliaReach/ReachabilityAnalysis.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/ReachabilityAnalysis.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/juliareach/ReachabilityAnalysis.jl?svg=true)](https://ci.appveyor.com/project/mforets/ReachabilityAnalysis-jl)
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

💾 [Installation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--installation)

📙 [Documentation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--documentation)

🎨 [Features](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--features)

🐾 [Examples Gallery](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--examples)

:blue_book: [Publications](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--blue_book-publications)

📜 [Citation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--citation)

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
| ![quadrotor](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-Quadrotor.png?raw=true) [Quadrotor altitude control](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Quadrotor/) |  ![LVHybrid](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-LotkaVolterra.png?raw=true) [Lotka-Volterra with tangential guard crossing](https://github.com/JuliaReach/ARCH2020_NLN_RE/blob/master/models/LotkaVolterra/lotka_volterra.jl)|
|         |      |
| ![LaubLoomis](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-LaubLoomis.png?raw=true) [Laub-Loomis model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/LaubLoomis/)    | ![PD](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-ProductionDestruction.png?raw=true)<br> [Production-Destruction model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/ProductionDestruction/)|
|         |      |
|![vanderpol](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-VanDerPol.png?raw=true) [Coupled van der pol oscillator model]() | ![spaccecraft](https://github.com/JuliaReach/JuliaReach-website/blob/master/src/images/NLN/ARCH-COMP20-JuliaReach-Spacecraft.png?raw=true) [Spacecraft rendez-vous model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Spacecraft/) |


## :blue_book: Publications

This library has been applied in a number of scientic works. We list them in reverse chronological order. 

[10] **Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.** Marcelo Forets, Daniel Freire, Christian Schilling, 2020. [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325).

[9] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Zongnan Bao, Marcelo Forets, Daniel Freire, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, Stefan Schupp, and Mark Wetzlinger (2020) ARCH20. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 16--48. [10.29007/7dt2](https://easychair.org/publications/paper/DRpS).

[8] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Pieter Collins, Marcelo Forets, Daniel Freire, Fabian Immler, Niklas Kochdumper, David P. Sanders and Christian
Schilling (2020) ARCH20. To appear in 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 49--75. [10.29007/zkf6](https://easychair.org/publications/paper/nrdD).

[7] **Case Study: Reachability Analysis of a unified Combat-Command-and-Control Model.** Sergiy Bogomolov, Marcelo Forets, Kostiantyn Potomkin. Accepted in [14th International Conference on Reachability Problems 2020](https://www.irif.fr/~rp2020/).

[6] **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. [Get pdf from arXiv: 1905.02458](https://arxiv.org/abs/1905.02458). Accepted in [Embedded Systems Week 2020](http://esweek.hosting2.acm.org/).

[5] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14--40 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[4] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41--61 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[3] **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Published in Proceedings of [HSCC'19](http://hscc2019.eecs.umich.edu/): 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'19), see [ACM link here](https://dl.acm.org/citation.cfm?id=3311804). [Get pdf from arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

[2] **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: 10.29007/73mb.

[1] **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) [HSCC'18](https://www.hscc2018.deib.polimi.it/) Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the [ACM Digital Library link](http://dx.doi.org/10.1145/3178126.3178128), or the [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526).

*Note:* Articles [1-7] use the former codebase `Reachability.jl`.

## 📜  Citation

Research credit and full references to the scientific papers presenting the algorithms implemented in this package can be found in the source code for each algorithm and in the [References](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/references/) section of the online documentation.

If you use this package for your research, we kindly ask you to cite the following paper, see [CITATION.bib](http://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/CITATION.bib). Moreover, please **also cite the appropriate original references to the algorithms used.**



