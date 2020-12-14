# ======================
# Dependencies
# ======================

using LinearAlgebra, SparseArrays, # modules from the Julia standard library
      Reexport,                    # see @reexport macro below
      RecipesBase,                 # plotting
      Parameters,                  # structs with kwargs
      StaticArrays,                # statically sized arrays
      RecursiveArrayTools          # vector of arrays
      #ExponentialUtilities        # (optional) Krylov subspace approximations

# the reexport macro ensures that the names exported by the following libraries
# are made available after loading ReachabilityAnalysis
@reexport using HybridSystems,
                IntervalMatrices,
                LazySets,
                MathematicalSystems,
                TaylorIntegration

# required to avoid conflicts with MathematicalSystems
using LazySets: LinearMap, AffineMap, ResetMap

# required to avoid conflicts with IntervalMatrices
using LazySets: Interval, isdisjoint, radius, sample, ∅, dim

# LazySets internal functions frequently used
using LazySets.Arrays: projection_matrix, SingleEntryVector

# aliases for intervals
const IM = IntervalMatrices
import IntervalArithmetic
const IA = IntervalArithmetic
const TimeInterval = IA.Interval{Float64}
import TaylorModels
const TM = TaylorModels

# method extensions for Taylor model reach-sets
import TaylorModels: domain, remainder, polynomial, get_order

@inline function _isapprox(Δt::TimeInterval, Δs::TimeInterval)
    return (inf(Δt) ≈ inf(Δs)) && (sup(Δt) ≈ sup(Δs))
end

# aliases for set types
const CPA = CartesianProductArray

# aliases for system types (TODO: merge with definitions in normalization.jl)
#const ACS = AbstractContinuousSystem
#const ADS = AbstractDiscreteSystem
#const LCS = LinearContinuousSystem
#const CLCCS = ConstrainedLinearControlContinuousSystem

# convenience union for dispatch on structs that are admissible as initial sets or inputs
const AdmissibleSet = Union{LazySet, UnionSet, UnionSetArray, IA.Interval, IA.IntervalBox}

# method extensions
import LazySets: dim, overapproximate, project, Projection,
                 linear_map, LinearMap, _split, split!, set, array

import Base: ∈, convert

# ======================
# Useful constants
# ======================

@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline unitBox(m) = IntervalBox(IA.Interval(0.0, 1.0), m)
@inline symBox(n::Integer) = IntervalBox(symI, n)
const zeroI = IA.Interval(0.0) # TODO use number type
const oneI = IA.Interval(1.0)
const symI = IA.Interval(-1.0, 1.0)

# ======================
# Optional dependencies
# ======================

using Requires

# convenience macro to annotate that a package is required
# usage:
# function foo(...)
#   @require MyPackage
#   ... # functionality that requires MyPackage to be loaded
# end
macro requires(module_name)
    m = Meta.quot(Symbol(module_name))
    return esc(:(@assert isdefined(@__MODULE__, $m) "package `$($m)` is required " *
                    "for this function; do `using $($m)` and try again"))
end

function __init__()
    @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include("init_DifferentialEquations.jl")
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" include("init_ExponentialUtilities.jl")
    @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" include("init_ModelingToolkit.jl")
end
