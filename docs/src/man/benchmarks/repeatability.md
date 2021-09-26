```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

## Repeatability evaluations

Traditionally, re-creation of computational results in research work is a challenging
task because details of the implementation are unavoidably absent in the paper.
Some authors post their code and data to their websites, but there is little formal
incentive to do so and no easy way to determine whether others can actually use the result.
As a consequence, computational results often become non reproducible -- even by the
research group which originally produced them -- after just a few years.

More recently, scientific conferences encourage and sometimes require that
authors improve the reproducibility of their computational results by providing suitable
*reproducibility evaluation* (RE) packages. Such RE are self-contained codes that
can be used to run and reproduce the results from the paper in a host machine.
Fortunately, the tools available for the Julia language are very convenenient
to prepare and distribute RE packages. The following RE packages are available:

- [ARCH2020 NLN RE](https://github.com/JuliaReach/ARCH2020_NLN_RE) -- Repeatability Evaluation package for the ARCH2020 Nonlinear Continuous and Hybrid Systems Competition.

- [ARCH2020 AFF RE](https://github.com/JuliaReach/ARCH2020_AFF_RE) -- Repeatability Evaluation package for the ARCH2020 Linear and Hybrid Systems Competition.

- [HSCC2019 RE](https://github.com/JuliaReach/HSCC2019_RE) -- Repeatability Evaluation (RE) package for the paper *JuliaReach: a Toolbox for Set-Based Reachability* published at the HSCC'2019 conference.

- [ARCH2019 RE](https://github.com/JuliaReach/ARCH2019_RE) -- Repeatability Evaluation package for the ARCH2019 Competition.

- [ARCH2018 RE](https://github.com/JuliaReach/ARCH2018_RE) -- Repeatability Evaluation package for the ARCH2018 Competition.

For installation instructions, see the `README` file in each package.
We have examples using virtual machines, and also using Docker containers for RE packages.
The advantage of using a Docker container is that downloading the necessary requirements
(including Julia itself) is an automated process.
