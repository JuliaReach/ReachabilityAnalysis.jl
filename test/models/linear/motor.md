```@meta
EditURL = "<unknown>/motor.jl"
```

## Model description

TODO: Add source for this model

The "motor" model is a linear system with input. It has 8 state variables
and 2 input variables.

```@example motor
using SparseArrays

I = [1, 2, 2, 3, 3, 3, 3, 4, 5, 6, 6, 7, 7, 7, 7, 8]
J = [2, 3, 2, 1, 2, 3, 4, 1, 6, 7, 6, 5, 6, 7, 8, 5]
vals = [1, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0, 1.0,
        8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0]
A = sparse(I, J, vals)
```

## Reachability settings

input set

```@example motor
B = sparse([4, 8], [1, 2], [-1.0, -1.0])
U = Hyperrectangle([0.23, 0.3], [0.07, 0.1])
```

Next, we instantiate the continuous LTI system,

```@example motor
S = @system(x' = Ax + Bu, u ∈ U, x ∈ Universe(8))

S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U);
nothing #hide
```

## Results

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

