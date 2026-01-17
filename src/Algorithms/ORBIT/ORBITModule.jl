module ORBITModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, ReachSet,
                              TimeInterval, compute_nsteps, homogenize, zeroI,
                              _normalize
using ..DiscretizationModule: next_set
using ..Exponentiation: _exp
using ..NoBloatingModule: NoBloating

using IntervalArithmetic: (..)
using LazySets: AbstractSingleton, LazySet, Singleton, Universe, ZeroSet, element, sample
using LinearAlgebra: mul!
using MathematicalSystems: AbstractContinuousSystem, AbstractInput, IVP,
                           discretize, initial_state, inputset, state_matrix,
                           stateset
using Parameters: @unpack
using ReachabilityBase.Arrays: vector_type
using Requires: @require

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("ORBIT.jl")
include("post.jl")

include("init.jl")

export ORBIT

end  # module
