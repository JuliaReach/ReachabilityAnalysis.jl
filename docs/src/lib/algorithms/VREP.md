```@docs
VREP
```

## Method

This algorithm solves the recurrence relation

```math
X_{k+1} = \Phi X_k \oplus V_k
```
using the vertex representation of the sets involved. If the system is homogeneous,
i.e. $V_k = ∅$ for all $k \geq 0$, then the number of vertices of the sequence doesn't
increase. On the other hand, the number of vertices of $X_k$ increases hence the method
requires doing some conservative reduction strategy.

The method is aimed towards small dimensional systems -- typically $n < 15$ --
which show a substantial benefit by using heap-allocated statically sized arrays
through the [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package.
To use static arrays, pass the option `static=true`. Optionally, a dimension field
can be passed to the `VREP` algorithm constructor too using the `dim` keyword argument.
(It is sometimes useful that Julia is able to infer the set representation, or type,
of the output flowpipe based only on the algorithm choice, but since the length of
a vector is stored as type information for the case of static arrays, it has to be given in advance.)

## Specifying the backend

If the dimension of the system is two, `VREP` uses efficient algorithms for convex
polygons (e.g. convex hull, Minkowski sum) implemented in
[LazySets.jl](https://github.com/JuliaReach/LazySets.jl/).
On the other hand, for systems of dimension higher than two, concrete polyhedral
computations use the [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)
libary which itself relies on specific *backends*, which
can be specified with the `backend` keyword argument in the `VREP` algorithm constructor.
Such backend is used in the discretization phase. Actually, `Polyhedra.jl` features a default
solver (hence it doesn't require additional packages apart from `Polyhedra.jl` itself),
although for efficiency reasons we suggest to use
[CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl).
For example, to use `CDDLib` in a 4-dimensional system, do

```julia
using ReachabilityAnalysis, Polyhedra, CDDLib

# ... define model ...
# prob = @ivp(...)
alg = VREP(δ=1e-3, static=true, dim=4, backend=CDDLib.Library()))
solve(prob, tspan=(0.0, 1.0), alg=alg)
```
Other backends are available e.g. [QHull.jl](https://github.com/JuliaPolyhedra/QHull.jl/).
