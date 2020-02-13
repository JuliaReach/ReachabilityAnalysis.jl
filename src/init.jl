# modules from the Julia standard library
using LinearAlgebra, SparseArrays

using Reexport, RecipesBase, Parameters

# the reexport macro ensures that the following libraries are made available upon
# loading this package
@reexport using LazySets,
                MathematicalSystems,
                HybridSystems,
                TaylorIntegration,
                IntervalMatrices

# required to avoid conflicts with MathematicalSystems
using LazySets: LinearMap, AffineMap, ResetMap, Interval

# required to avoid conflicts with IntervalMatrices and IntervalArithmetic
using LazySets: isdisjoint, radius, sample, ∅, dim

# aliases for intervals
const IM = IntervalMatrices
import IntervalArithmetic
const IA = IntervalArithmetic
const I64 = IA.Interval{Float64}

# aliases for sets / system types
const CPA = CartesianProductArray
const LCS = LinearContinuousSystem
const CLCCS = ConstrainedLinearControlContinuousSystem
