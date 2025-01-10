# Systems

This section of the manual describes the systems types that are used in this library.

[MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl)
provides some convenience types and methods to work with mathematical systems models. Every system inherits
from `AbstractSystem`.

## Linear systems

Two commonly used types of systems are discrete and continuous systems.

**Discrete systems.** A discrete system consists of a matrix representing the system dynamics, a set
of initial states, a set of nondeterministic inputs, and a discretization step
δ.

**Continuous systems.**  A continuous system consists of a matrix representing the system dynamics, a set
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

## Second order systems

A second order system is one of the form

```math
    Mx''(t) + Cx'(t) + Kx(t) = f(t)
```
where ``x(t) ∈ \mathbb{R}^n`` is the state vector and ``f : \mathbb{R} \to \mathbb{R}^n``
is the forcing term. Here ``M``, ``C`` and ``K`` are often called the mass matrix,
viscosity matrix and stiffness matrix respectively. These names are adopted from
physical applications, particularly from structural mechanics. Assuming that the
matrix ``M`` is invertible, we can transform the second order system to a first order
system introducing auxiliary variables, ``x̃(t) = [x(t),~v(t)]^T``, where ``v(t) := x'(t)`` is the vector of velocities. Then,

```math
    x̃(t)' = Ax̃(t) + Bf(t)
```
where

```math
A = \begin{pmatrix}
0 && I \\ -M^{-1}K && -M^{-1}C
\end{pmatrix},\qquad B = \begin{pmatrix}
0  \\ M^{-1}
\end{pmatrix}
```
See the [SecondOrder](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.SecondOrderConstrainedLinearControlContinuousSystem) documentation
in `MathematicalSystems.jl` for additional details in second order ODEs types.

!!! note
    A similar relation can be obtained using the alternative convention ``[v(t),~x(t)]^T``.
    Use `derivatives_first=true` in the `normalize` function to swap between these conventions
    (it is set to `false` by default).

## Normalization

## Homogenization

## Nonlinear systems

## Parametric systems
