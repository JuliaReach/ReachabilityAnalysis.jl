```@docs
VREP
```

## Specifing the backend

If the dimension of the system is two, `VREP` uses efficient algorithms for convex
polygons (e.g. convex hull, Minkowski sum) implemented in `LazySets.jl`. On the other
hand, for systems of dimension higher than two the polyhedral computations backend
can be specified with the `backend` keyword argument. To do so, load `Polyhedra`
and optionally a backend such as `CDDLib` or `QHull` and let
`VREP(δ=1e-3, static=true, dim=4, backend=CDDLib.Library()))`.
