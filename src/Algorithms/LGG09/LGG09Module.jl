module LGG09Module

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe,
                              TemplateReachSet, TimeInterval, VecOrTuple,
                              hasinput, zeroI, _collect, _isdisjoint
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
using Requires: @require
using SparseArrays: sparse

import ..ReachabilityAnalysis: post, numtype, rsetrep, setrep, step_size

include("LGG09.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

include("init.jl")

export LGG09

end  # module
