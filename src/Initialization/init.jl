# ======================
# Dependencies
# ======================

using LinearAlgebra, SparseArrays, # modules from the Julia standard library
      Reexport,                    # see @reexport macro below
      RecipesBase,                 # plotting
      Parameters,                  # structs with kwargs
      StaticArrays,                # statically sized arrays
      RecursiveArrayTools          # vector of arrays type

# manipulate function definition expressions
using ExprTools: splitdef, combinedef

# the reexport macro ensures that the names exported by the following libraries
# are made available after loading ReachabilityAnalysis
@reexport using HybridSystems,
                IntervalMatrices,
                ReachabilityBase,
                LazySets,
                MathematicalSystems,
                TaylorIntegration

# required to avoid conflicts with MathematicalSystems
using LazySets: AffineMap, ResetMap

# required to avoid conflicts with IntervalMatrices
using LazySets: Interval, radius, sample, ∅, dim, scale, scale!

# in-place set operations
using LazySets: linear_map!

# LazySets internal functions frequently used
using ReachabilityBase.Arrays: projection_matrix, SingleEntryVector,
                               isinvertible, vector_type
using LazySets.Approximations: AbstractDirections
using LazySets: @commutative, AbstractReductionMethod

# aliases for intervals
const IM = IntervalMatrices
import IntervalArithmetic
const IA = IntervalArithmetic
const TimeInterval = IA.Interval{Float64}
import TaylorModels
const TM = TaylorModels

# method extensions for Taylor model reach-sets
import TaylorModels: domain, remainder, polynomial, get_order, evaluate

# aliases for set types
const CPA = CartesianProductArray

# convenience union for dispatch on structs that are admissible as initial sets or inputs
const AdmissibleSet = Union{LazySet,UnionSet,UnionSetArray,IA.Interval,IA.IntervalBox}

# method extensions
import LazySets: dim, overapproximate, box_approximation, project, Projection,
                 intersection, directions, linear_map, LinearMap, _split, split!,
                 set, array, _isapprox

import MathematicalSystems: discretize
import Base: ∈, ∩, convert, isdisjoint
import LinearAlgebra: normalize

import CommonSolve: solve # common solve name

# ======================
# Useful constants
# ======================

@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline unitBox(m) = IntervalBox(IA.interval(0.0, 1.0), m)
@inline symBox(n::Integer) = IntervalBox(symI, n)
const zeroI = IA.interval(0.0) # TODO use number type
const oneI = IA.interval(1.0)
const symI = IA.interval(-1.0, 1.0)

# common aliases for system's names
const LCS = LinearContinuousSystem
const LDS = LinearDiscreteSystem
const CLCS = ConstrainedLinearContinuousSystem
const CLDS = ConstrainedLinearDiscreteSystem
const CLCCS = ConstrainedLinearControlContinuousSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem
const ACS = AffineContinuousSystem
const ADS = AffineDiscreteSystem
const CACCS = ConstrainedAffineControlContinuousSystem
const CACDS = ConstrainedAffineControlDiscreteSystem
const CACS = ConstrainedAffineContinuousSystem
const CADS = ConstrainedAffineDiscreteSystem
const BBCS = BlackBoxContinuousSystem
const CBBCS = ConstrainedBlackBoxContinuousSystem
const CBBCCS = ConstrainedBlackBoxControlContinuousSystem
const SOLCS = SecondOrderLinearContinuousSystem
const SOACS = SecondOrderAffineContinuousSystem
const SOCLCCS = SecondOrderConstrainedLinearControlContinuousSystem
const SOCACCS = SecondOrderConstrainedAffineControlContinuousSystem
const SecondOrderSystem = Union{SOLCS,SOACS,SOCLCCS,SOCACCS}
const NonlinearSystem = Union{BBCS,CBBCS,CBBCCS}

@inline function _isapprox(Δt::TimeInterval, Δs::TimeInterval)
    return (inf(Δt) ≈ inf(Δs)) && (sup(Δt) ≈ sup(Δs))
end

const VecOrTuple = Union{<:AbstractVector{Int},NTuple{D,Int}} where {D}
const VecOrTupleOrInt = Union{<:AbstractVector{Int},NTuple{D,Int},Int} where {D}

# ======================
# Optional dependencies
# ======================

using Requires

function __init__()
    # numerical differential equations suite
    @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include("init_DifferentialEquations.jl")

    # sparse dynamic representation of multivariate polynomials
    @require DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07" include("init_DynamicPolynomials.jl")

    # exponentiation methods using Krylov subspace and Padé approximations
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" include("init_Expokit.jl")

    # exponentiation methods using Krylov subspace approximations
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" include("init_ExponentialUtilities.jl")

    # exponentiation methods for large sparse matrices
    @require FastExpm = "7868e603-8603-432e-a1a1-694bd70b01f2" include("init_FastExpm.jl")

    # external reachability backend: Flow* C++ library
    @require Flowstar = "a8054ddd-9dca-4d20-8ffe-ae96ec1541f1" include("init_Flowstar.jl")

    # methods that require optimization modeling and solvers
    @require JuMP = "4076af6c-e467-56ae-b986-b466b2749572" include("init_JuMP.jl")

    # tools for symbolic algebra
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")

    # tools for symbolic computation
    #@require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" include("init_ModelingToolkit.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include("init_Symbolics.jl")
end

# ===========================
# Utility macros
# ===========================

"""
    @requires(module_name)

Convenience macro to annotate that a package is required to use a certain function.

### Input

- `module_name` -- name of the required package

### Output

The macro expands to an assertion that checks whether the module `module_name` is
known in the calling scope.

### Notes

Usage:

```julia
function foo(...)
    @require MyPackage
    ... # functionality that requires MyPackage to be loaded
end
```
"""
macro requires(module_name)
    m = Meta.quot(Symbol(module_name))
    return esc(:(@assert isdefined(@__MODULE__, $m) "package `$($m)` is required " *
                                                    "for this function; do `using $($m)` and try again"))
end
