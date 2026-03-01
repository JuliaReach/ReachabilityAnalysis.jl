module BFFPSV18Module

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, SparseReachSet,
                              TimeInterval, hasinput, zeroT, _isdisjoint
using ..ForwardModule: Forward

using IntervalArithmetic.Symbols: (..)
using LazySets: CartesianProductArray, Hyperrectangle, Interval, LazySet,
                LinearMap, MinkowskiSumArray, Projection, Universe, decompose,
                overapproximate, ⊕
using LinearAlgebra: mul!
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, stateset
using Parameters: @unpack
using SparseArrays: SparseMatrixCSC, sparse

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("BFFPSV18.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export BFFPSV18

end  # module
