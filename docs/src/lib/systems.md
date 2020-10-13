## Representation of linear systems

[MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl)
provides some convenience types and methods to work with mathematical systems models. Every system inherits
from `AbstractSystem`.

We support the following two concrete types of systems.

### Discrete system

A discrete system consists of a matrix representing the system dynamics, a set
of initial states, a set of nondeterministic inputs, and a discretization step
δ.

### Continuous system

A continuous system consists of a matrix representing the system dynamics, a set
of initial states, and a set of nondeterministic inputs.

## Nondeterministic inputs

The above systems may contain nondeterministic inputs, which are wrapped in
special types. Every nondeterministic input representation inherits from
`NonDeterministicInput`.

The inputs are closely related to a `DiscreteSystem` in the sense that
for each discrete time step the input set may change. We support iteration
through the inputs over time.

### Constant nondeterministic inputs

Constant nondeterministic inputs are chosen from a set of values that does not
change over time. Note that, while the set is constant, the inputs themselves
vary over time.

### Time-varying nondeterministic inputs

Time-varying nondeterministic inputs are chosen from a set of values that
changes over time (with each time step).

## Normalization
