# About

This page contains some general information about this project, and recommendations
about contributing.

```@contents
Pages = ["about.md"]
```

## Contributing

If you like this package, consider contributing! We welcome bug reports,
examples that can be added to the documentation, or new feature proposals.

Below some conventions that we follow when contributing
to this package are detailed. For specific guidelines on documentation, see the [Documentations Guidelines wiki](https://github.com/JuliaReach/LazySets.jl/wiki/Documentation-Guidelines).

#### Branches

Each pull request (PR) should be pushed in a new branch with the name of the author
followed by a descriptive name, e.g. `t/mforets/my_feature`. If the branch is associated
to a previous discussion in one issue, we use the name of the issue for easier
lookup, e.g. `t/mforets/7`.

### Unit testing and continuous integration (CI)

This project is synchronized with Travis CI, such that each PR gets tested
before merging (and the build is automatically triggered after each new commit).
For the maintainability of this project, it is important to understand and fix the
failing doctests if they exist. We develop in Julia v0.6.0, but for experimentation
we also build on the nightly branch.

To run the unit tests locally, you should do:

```julia
$ julia --color=yes test/runtests.jl
```

### Contributing to the documentation

This documentation is written in Markdown, and it relies on
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to produce the HTML
layout. To build the docs, run `make.jl`:

```julia
$ julia --color=yes docs/make.jl
```

## References

This repository was originally motivated by the mathematical approach described in *Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices*,  Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Andreas Podelski, Christian Schilling, Frédéric Viry, in [21st ACM International Conference on Hybrid Systems: Computation and Control](https://www.hscc2018.deib.polimi.it/), 2018 Edition (Porto, Portugal), see the [arXiv pre-print here](https://arxiv.org/abs/1801.09526).

For a full references list of algorithms implemented in this repository, consult [the References wiki](https://github.com/JuliaReach/Reachability.jl/wiki/References).


## Developers

Here we list the names of the contributors to `ReachabilityAnalysis.jl`
(in alphabetic order).

- [Marcelo Forets](http://github.com/mforets)
- [Daniel Freire](http://github.com/dfcaporale)
- [Christian Schilling]http://github.com/schillic)

### Acknowledgements

We are also grateful to the following persons for enlightening discussions:

- [Sergiy Bogomolov](https://www.sergiybogomolov.com/)
- [Goran Frehse](https://sites.google.com/site/frehseg/)
- Nikolaos Kekatos
- Andreas Podelski
- Kostiantyn Potomkin
- Alexandre Rocca
- Frederic Viry
