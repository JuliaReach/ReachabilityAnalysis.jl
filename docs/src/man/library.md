# Model library

The library [ReachabilityModels.jl](https://github.com/JuliaReach/ReachabilityModels.jl)
contains a collection of pre-made models that can be found in books, articles or
other publicly available software related to reachability analysis. As it is
explained in the [documentation of that library](https://juliareach.github.io/ReachabilityModels.jl/dev/),
once installed use it as:

```julia
using ReachabilityModels, Plots

prob = fetch_model("building") # initial-value problem

sol = solve(prob, T=5.0); # solve it using default options

plot(sol, vars=(0, 25)) # plot the solution
```
