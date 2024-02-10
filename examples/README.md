This folder contains the source code for the examples available in the
[Examples](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/examples) section of the
[online documentation](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/).
The examples are also available in the form of Jupyter notebooks, and can be
visualized online with [nbviewer](https://nbviewer.jupyter.org/).
The examples are generated using the package [Literate.jl](https://github.com/fredrikekre/Literate.jl).

See also the library [ReachabilityModels.jl](https://github.com/JuliaReach/ReachabilityModels.jl),
which contains additional example models, and implements a convenient API
to load models that can be used like this:

```julia
using ReachabilityModels, Plots

prob = load_model("building"); # initial-value problem

sol = solve(prob, tspan=(0, 5)); # solve it using default options

plot(sol, vars=(0, 25)) # plot the projected flowpipe of x25 vs time
```
