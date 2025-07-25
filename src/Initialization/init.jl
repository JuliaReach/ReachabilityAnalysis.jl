# ======================
# Dependencies
# ======================

using LinearAlgebra, SparseArrays, # modules from the Julia standard library
      Reexport,                    # see @reexport macro below
      RecipesBase,                 # plotting
      Parameters,                  # structs with kwargs
      StaticArrays                 # statically sized arrays

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
using LazySets: Interval, radius, sample, ∅, dim, scale, scale!, ⊂, matrix

# JuliaReach internal functions
using ReachabilityBase.Arrays: projection_matrix, SingleEntryVector,
                               isinvertible, samedir, vector_type
using ReachabilityBase.Comparison: _leq, _geq
using ReachabilityBase.Require: @required
using LazySets.Approximations: AbstractDirections
using LazySets: @commutative, AbstractReductionMethod, linear_map!

# aliases for intervals
const IM = IntervalMatrices
import IntervalArithmetic as IA
const TimeInterval = IA.Interval{Float64}
import TaylorModels as TM
using TaylorModels: TaylorModel1, TaylorN, fp_rpa, shrink_wrapping!

# method extensions for Taylor model reach-sets
import TaylorModels: domain, remainder, polynomial, get_order, evaluate

# aliases for set types
const CPA = CartesianProductArray

# convenience union for dispatch on structs that are admissible as initial sets or inputs
const AdmissibleSet = Union{LazySet,IA.Interval,IA.IntervalBox}

# method extensions
import LazySets: dim, overapproximate, box_approximation, project, Projection,
                 intersection, directions, linear_map, LinearMap, split!,
                 set, array, _isapprox, flatten,
                 _plot_singleton_list_1D, _plot_singleton_list_2D

import Base: ∈, ∩, convert, isdisjoint
import LinearAlgebra: normalize

import MathematicalSystems: system, statedim, initial_state
import HybridSystems: HybridSystem

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
    @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" include("init_OrdinaryDiffEq.jl")

    # sparse dynamic representation of multivariate polynomials
    @require DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07" include("init_DynamicPolynomials.jl")

    # exponentiation methods using Krylov subspace approximations
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" include("init_ExponentialUtilities.jl")

    # external reachability backend: Flow* C++ library
    @require Flowstar = "a8054ddd-9dca-4d20-8ffe-ae96ec1541f1" include("init_Flowstar.jl")

    # tools for symbolic algebra
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")

    # tools for symbolic computation
    #@require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" include("init_ModelingToolkit.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include("init_Symbolics.jl")
end
