# About

This page contains some general information about this project, and recommendations about contributing.


## Contributing

If you like this package, consider contributing! We welcome bug reports,
examples that can be added to the documentation, or new feature proposals.

Below we list some conventions that we follow when contributing
to this package. For specific guidelines on documentation see the
[JuliaReach Developer's Documentation](https://github.com/JuliaReach/JuliaReachDevDocs).

#### Branches

Each pull request (PR) should be pushed in a new branch with the name of the author
followed by a descriptive name, e.g. `mforets/my_feature`. If the branch is associated
to a previous discussion in an issue, we use the name of the issue for easier
lookup, e.g. `mforets/7`.

### Unit testing and continuous integration (CI)

This project is synchronized with Travis CI, such that each PR gets tested
before merging (and the build is automatically triggered after each new commit).
For the maintainability of this project, it is important to understand and fix the
failing tests if they exist. We support Julia v1.0 and develop in the latest stable release.
For experimentation we also build on the nightly branch.

To run the unit tests locally, you should do:

```julia
$ julia --color=yes test/runtests.jl
```

### Contributing to the documentation

This documentation is written in Markdown, and it relies on
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to produce the HTML
layout. To build the documentation, run `make.jl`:

```julia
$ julia --color=yes docs/make.jl
```
