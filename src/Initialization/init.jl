# ======================
# Dependencies
# ======================

using LinearAlgebra, SparseArrays, # modules from the Julia standard library
      Reexport,    # see @reexport macro below
      RecipesBase, # plotting
      Parameters,   # structs with kwargs
      StaticArrays,
      RecursiveArrayTools

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
end
