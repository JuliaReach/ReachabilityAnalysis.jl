# ======================
# Dependencies
# ======================

using LinearAlgebra: Diagonal, I, dot, isdiag, mul!, tr
using Parameters: @unpack, @with_kw
using RecipesBase: @recipe, @series
using Reexport: @reexport
using SparseArrays: SparseMatrixCSC, SparseVector, sparsevec, spzeros
using StaticArrays: MMatrix, SMatrix, SVector

# the reexport macro ensures that the names exported by the following libraries
# are made available after loading ReachabilityAnalysis
@reexport using HybridSystems,
                IntervalMatrices,
                LazySets,
                MathematicalSystems,
                TaylorIntegration

# required to avoid conflicts with MathematicalSystems
using LazySets: AffineMap, LinearMap, ResetMap

# required to avoid conflicts with IntervalMatrices
using LazySets: Interval, radius, sample, ∅, dim, scale, scale!, ⊂, matrix

# JuliaReach internal functions
import ReachabilityBase
using ReachabilityBase.Arrays: projection_matrix, SingleEntryVector,
                               isinvertible, samedir
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: _leq, _geq
using ReachabilityBase.Require: @required
using LazySets.Approximations: AbstractDirections
using LazySets: AbstractReductionMethod, linear_map!

# aliases for intervals
const IM = IntervalMatrices
import IntervalArithmetic as IA
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
                 intersection, directions, linear_map, split!, set, array,
                 constrained_dimensions
import ReachabilityBase.Comparison: _isapprox
import Base: +, -, *, ∈, ∩, ⊆, convert, isdisjoint, isempty
import LinearAlgebra: normalize

import MathematicalSystems: system, statedim, initial_state
import HybridSystems: HybridSystem, guard, mode

import CommonSolve: solve # common solve name

# ======================
# Useful constants
# ======================

@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline symBox(n::Integer) = IntervalBox(symI, n)
const zeroI = IA.interval(0.0) # TODO use number type
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
const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem

const VecOrTuple = Union{<:AbstractVector{Int},NTuple{D,Int}} where {D}
const VecOrTupleOrInt = Union{<:AbstractVector{Int},NTuple{D,Int},Int} where {D}

# ======================
# Optional dependencies
# ======================

using Requires: @require

function __init__()
    # numerical differential equations suite
    @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" include("init_OrdinaryDiffEq.jl")

    # sparse dynamic representation of multivariate polynomials
    @require DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07" include("init_DynamicPolynomials.jl")

    # tools for symbolic algebra
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")

    # tools for symbolic computation
    #@require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" include("init_ModelingToolkit.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include("init_Symbolics.jl")
end
