"""
Interface for conservative time discretization methods.
"""
module DiscretizationModule

using MathematicalSystems
import IntervalArithmetic as IA
using IntervalMatrices
using LazySets
using LazySets.Approximations: AbstractDirections
using Reexport

using ..Exponentiation
import ..Exponentiation: _alias

@reexport import MathematicalSystems: discretize
export AbstractApproximationModel, sih, isinterval, next_set

# convenience functions
next_set(inputs::ConstantInput) = collect(nextinput(inputs, 1))[1]
next_set(inputs::AbstractInput, state::Int64) = collect(nextinput(inputs, state))[1]

abstract type AbstractApproximationModel end

function _default_approximation_model(::IVP{<:AbstractContinuousSystem})
    return Forward()
end

# some algorithms require a polyhedral computations backend
hasbackend(alg::AbstractApproximationModel) = false

# symmetric interval hull options
sih(X, ::Val{:lazy}) = SymmetricIntervalHull(X)
sih(X, ::Val{:concrete}) = LazySets.symmetric_interval_hull(X)

# interval matrix functions
isinterval(A::AbstractMatrix{N}) where {N<:Number} = false
isinterval(A::IntervalMatrix{N,IT}) where {N,IT<:IA.Interval{N}} = true
isinterval(A::AbstractMatrix{IT}) where {IT<:IA.Interval} = true

# options for a-posteriori transformation of a discretized set
_alias(setops::AbstractDirections) = setops
_alias(setops::Val{:lazy}) = setops
_alias(setops::Val{:concrete}) = setops
_alias(setops::Val{:vrep}) = setops
_alias(setops::Val{:box}) = setops
_alias(setops::Val{:zono}) = setops
_alias(::Val{:zonotope}) = Val(:zono)

"""
    discretize(ivp::IVP, δ, alg::AbstractApproximationModel)

Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.

### Input

- `ivp`   -- initial value problem for a linear ODE in canonical form (see `Notes` below)
- `δ`     -- step size
- `alg`   -- algorithm used to compute the approximation model

### Output

The initial value problem of a discrete system.

### Notes

Different approximation algorithms and their respective options are described
in the docstring of each method. Here is a list of all the available approximation models:

```jldoctest
julia> subtypes(ReachabilityAnalysis.DiscretizationModule.AbstractApproximationModel)
9-element Vector{Any}:
 Backward
 CorrectionHull
 CorrectionHullMatrixZonotope
 FirstOrder
 FirstOrderZonotope
 Forward
 ForwardBackward
 NoBloating
 SecondOrderddt
 StepIntersect
```

Initial-value problems considered in this function are of the form

```math
x' = Ax(t) + u(t),\\qquad x(0) ∈ \\mathcal{X}_0,\\qquad (1)
```
and where ``u(t) ∈ U(k)`` add where ``\\{U(k)\\}_k`` is a sequence of sets of
non-deterministic inputs and ``\\mathcal{X}_0`` is the set of initial
states. Other problems, e.g. ``x' = Ax(t) + Bu(t)`` can be brought
to the canonical form with the function [`normalize`](@ref).

For references to the original papers introducing each algorithm, see the docstrings,
e.g. `?Forward`.
"""
function discretize(ivp::IVP, δ, alg::AbstractApproximationModel)
    return error("discretization not implemented for the given arguments: $ivp, $alg")
end

end  # module
