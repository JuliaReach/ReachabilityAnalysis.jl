module LGG09Module

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe,
                              TemplateReachSet, TimeInterval, VecOrTuple,
                              hasinput, zeroI, _is_intersection_empty
using ..ForwardModule: Forward

using Base: OneTo, Slice
using IntervalArithmetic: (..)
using LazySets: LazySet, Universe, dim, set, ρ
using LazySets.Approximations: AbstractDirections, BoxDirections,
                               CustomDirections, OctDirections
using LinearAlgebra: mul!
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, statedim, stateset
using Parameters: @unpack
using ReachabilityBase.Arrays: SingleEntryVector
using SparseArrays: sparse

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("LGG09.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export LGG09

end  # module
