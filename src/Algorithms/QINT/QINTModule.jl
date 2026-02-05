module QINTModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, HybridFlowpipe,
                              INT, ReachSet, StateInLocation, TimeInterval,
                              WaitingList, hasinput, state, tend, tspan, tstart,
                              zeroT, _normalize, _TimeInterval
using ..ForwardModule: Forward

using IntervalArithmetic: mid
using LazySets: ConvexHullArray, Interval, overapproximate, set, ⊕
using MathematicalSystems: AbstractContinuousSystem, AffineContinuousSystem,
                           IVP, initial_state, statedim
using Parameters: @unpack

import ..ReachabilityAnalysis: post, numtype, rsetrep, setrep, step_size

include("QINT.jl")
include("post.jl")
include("reach_homog.jl")
#include("reach_inhomog.jl")

export QINT

end  # module
