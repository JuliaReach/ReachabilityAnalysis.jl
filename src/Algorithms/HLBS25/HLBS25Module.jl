module HLBS25Module

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, IDGenerator,
                              ReachSet, TimeInterval, compute_nsteps, discretize,
                              fresh!, hasinput, synchronize!, tend, zeroI, zeroT
using ..CorrectionHullMatrixZonotopeModule: CorrectionHullMatrixZonotope,
                                            overapproximate_discrete_input_split

using IntervalArithmetic.Symbols: (..)
using LazySets: AbstractReductionMethod, AbstractZonotope, ExactSum,
                ExponentialMap, GIR05, MatrixZonotope, MatrixZonotopeExp,
                SparsePolynomialZonotope, Zonotope, center, concretize, dim,
                exact_sum, expmat, genmat, genmat_dep, genmat_indep, indexvector,
                minkowski_sum, overapproximate, reduce_order,
                remove_redundant_generators, scale
using LinearAlgebra: norm
using MathematicalSystems: AbstractContinuousSystem, AbstractDiscreteSystem, IVP,
                           initial_state, input_matrix, inputset, state_matrix
using Parameters: @unpack

import ..ReachabilityAnalysis: post, numtype, rsetrep, setrep, step_size

include("HLBS25.jl")
include("post.jl")
include("common.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export HLBS25

end  # module
