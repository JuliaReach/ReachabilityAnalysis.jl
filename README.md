# ReachabilityAnalysis.jl

[![Build Status](https://github.com/JuliaReach/ReachabilityAnalysis.jl/workflows/CI/badge.svg)](https://github.com/JuliaReach/ReachabilityAnalysis.jl/actions?query=workflow%3ACI)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/juliareach/ReachabilityAnalysis.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/juliareach/ReachabilityAnalysis.jl/coverage.svg?branch=master)](https://codecov.io/github/juliareach/ReachabilityAnalysis.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## ✨  What is Reachability Analysis?

Reachability analysis is concerned with computing rigorous approximations of the set
of states reachable by a dynamical system. In the scope of this package are systems
modeled by **continuous** or **hybrid** dynamical systems, where the dynamics changes with discrete events.
Systems are modelled by ordinary differential equations (ODEs) or semi-discrete partial differential equations (PDEs),
with uncertain initial states, uncertain parameters or non-deterministic inputs.

The library is oriented towards a class of numerical methods known as **set propagation techniques:**
to compute the set of states reachable by continuous or hybrid systems, such methods iteratively
propagate a sequence of sets starting from the set of initial states, according to the systems' dynamics.

See our [Commonly asked questions](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/man/faq/#Commonly-asked-questions-(FAQ)-1) (FAQ) section for pointers to the relevant literature, related tools and more.

## 🎯  Table of Contents

💾 [Installation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--installation)

📙 [Documentation](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--documentation) << WORK-IN-PROGRESS

🎨 [Features](https://github.com/JuliaReach/ReachabilityAnalysis.jl#--features)

:checkered_flag: [Quickstart](https://github.com/JuliaReach/ReachabilityAnalysis.jl#checkered_flag-quickstart)

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

The following types of systems are supported (click on the left arrow to display a list of examples):

<details>
  <summary> Continuous ODEs with linear dynamics :heavy_check_mark: </summary>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/OpAmp/">Operational amplifier</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/heat/">Heat</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/ISS/">ISS</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/motor">Motor</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Building/">Building</a> </p>
</details>


<details>
  <summary> Continuous ODEs with non-linear dynamics :heavy_check_mark: </summary>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Quadrotor/">Quadrotor</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Brusselator/">Brusselator</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/SEIR/">SEIR model</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/robot_arm">Robot arm</a> </p>
</details>

<details>
  <summary> Continuous ODEs with parametric uncertainty :heavy_check_mark: </summary>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/TransmissionLine/">Transmission line</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/LotkaVolterra/">Lotka-Volterra</a> </p>
</details>

<details>
  <summary> Hybrid systems with piecewise-affine dynamics :heavy_check_mark: </summary>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Platoon/">Platooning</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/bouncing_ball">Bouncing ball</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/navigation_system">Navigation system</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/thermostat">Thermostat</a> </p>
</details>

<details>
  <summary> Hybrid systems with non-linear dynamics :heavy_check_mark: </summary>
  <p> <a href="https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Spacecraft/">Spacecraft</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/cardiac_cell">Cardiatic cell</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/powertrain_control">Powetrain control</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/spiking_neuron">Spiking neuron</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/bouncing_ball_nonlinear">Bouncing ball</a> </p>
</details>

<details>
  <summary> Hybrid systems with clocked linear dynamics :heavy_check_mark: </summary>
  <p> <a href="https://github.com/JuliaReach/ARCH2020_AFF_RE/blob/master/models/EMBrake/embrake.jl">Electromechanic break</a> </p>
  <p> <a href="https://juliareach.github.io/ReachabilityModels.jl/dev/models/clocked_thermostat">Clocked thermostat</a> </p>
</details>

## :checkered_flag: Quickstart

In less than 15 lines of code, we can formulate, solve and visualize the set of states reachable by the [Duffing oscillator](https://en.wikipedia.org/wiki/Duffing_equation) starting from any initial condition
with position in the interval `0.9 .. 1.1` and velocity in `-0.1 .. 0.1`.

```julia
using ReachabilityAnalysis, Plots

const ω = 1.2
const T = 2*pi / ω

@taylorize function duffing!(du, u, p, t)
    local α = -1.0
    local β = 1.0
    local δ = 0.3
    local γ = 0.37
    
    x, v = u
    f = γ * cos(ω * t)

    # write the nonlinear differential equations defining the model
    du[1] = u[2]
    du[2] = -α*x - δ*v - β*x^3 + f
end

# set of initial states
X0 = Hyperrectangle(low=[0.9, -0.1], high=[1.1, 0.1])

# formulate the initial-value problem
prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim=2)

# solve using a Taylor model set representation
sol = solve(prob, tspan=(0.0, 20*T), alg=TMJets())

# plot the flowpipe in state-space
plot(sol, vars=(1, 2), xlab="x", ylab="v", lw=0.5, color=:red)
```

<img src="https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/duffing2.png?raw=true" alt="" width="600">


## 🐾  Examples

|         |      |
|:--------:|:-----:|
| ![quadrotor](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-Quadrotor.png?raw=true) [Quadrotor altitude control](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Quadrotor/) |  ![LVHybrid](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-LotkaVolterra.png?raw=true) [Lotka-Volterra with tangential guard crossing](https://github.com/JuliaReach/ARCH2020_NLN_RE/blob/master/models/LotkaVolterra/lotka_volterra.jl)|
|         |      |
| ![LaubLoomis](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-LaubLoomis.png?raw=true) [Laub-Loomis model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/LaubLoomis/)    | ![PD](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-ProductionDestruction.png?raw=true)<br> [Production-Destruction model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/ProductionDestruction/)|
|         |      |
|![vanderpol](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-VanDerPol.png?raw=true) [Coupled van der pol oscillator model](https://github.com/JuliaReach/ARCH2020_NLN_RE/blob/master/models/VanDerPol/vanderpol.jl) | ![spaccecraft](https://github.com/JuliaReach/JuliaReach-website/blob/master/images/reachability/NLN/ARCH-COMP20-JuliaReach-Spacecraft.png?raw=true) [Spacecraft rendez-vous model](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/models/Spacecraft/) |


## :blue_book: Publications

This library has been applied in a number of scientic works. We list them in reverse chronological order.

[11] **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems.** Forets, Marcelo, Daniel Freire Caporale, and Jorge M. Pérez Zerpa.. arXiv preprint [arXiv:2105.05841](https://arxiv.org/abs/2105.05841) (2021).

[10] **Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.** Marcelo Forets, Daniel Freire, Christian Schilling, 2020. [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325). In [18th ACM-IEEE International Conference on Formal Methods and Models for System Design
](https://iitjammu.ac.in/conferences/memocode2020/index.html).

[9] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Zongnan Bao, Marcelo Forets, Daniel Freire, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, Stefan Schupp, and Mark Wetzlinger (2020) ARCH20. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 16--48. [10.29007/7dt2](https://easychair.org/publications/paper/DRpS).

[8] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Pieter Collins, Marcelo Forets, Daniel Freire, Fabian Immler, Niklas Kochdumper, David P. Sanders and Christian
Schilling (2020) ARCH20. To appear in 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 49--75. [10.29007/zkf6](https://easychair.org/publications/paper/nrdD).

[7] **Case Study: Reachability Analysis of a unified Combat-Command-and-Control Model.** Sergiy Bogomolov, Marcelo Forets, Kostiantyn Potomkin. *International Conference on Reachability Problems (2020). Lecture Notes in Computer Science, vol 12448.* (2020) doi: [10.1007/978-3-030-61739-4_4](https://dx.doi.org/10.1007/978-3-030-61739-4_4). Presented in the [14th International Conference on Reachability Problems 2020](https://www.irif.fr/~rp2020/). [article](https://link.springer.com/chapter/10.1007/978-3-030-61739-4_4)

[6] **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. *IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 39:11 (2020).* doi: [10.1109/TCAD.2020.3012859](https://dx.doi.org/10.1109/TCAD.2020.3012859). Presented in [Embedded Systems Week 2020](http://esweek.hosting2.acm.org/). [Get pdf from arXiv: 1905.02458](https://arxiv.org/abs/1905.02458).

[5] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14--40 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[4] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41--61 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[3] **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Published in Proceedings of [HSCC'19](http://hscc2019.eecs.umich.edu/): 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'19), see [ACM link here](https://dl.acm.org/citation.cfm?id=3311804). [Get pdf from arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

[2] **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: [10.29007/73mb](https://dx.doi.org/10.29007/73mb).

[1] **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) [HSCC'18](https://www.hscc2018.deib.polimi.it/) Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the [ACM Digital Library link](http://dx.doi.org/10.1145/3178126.3178128), or the [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526).

*Note:* Articles [1-7] use the former codebase `Reachability.jl`.

## 📜  Citation

Research credit and full references to the scientific papers presenting the algorithms implemented in this package can be found in the source code for each algorithm and in the [References](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/references/) section of the online documentation.

If you use this package for your research, we kindly ask you to cite the following paper, see [CITATION.bib](http://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/CITATION.bib). Moreover, please **also cite the appropriate original references to the algorithms used.**
