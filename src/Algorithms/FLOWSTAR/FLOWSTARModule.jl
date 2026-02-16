module FLOWSTARModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, TimeInterval,
                              TaylorModelReachSet, domain, tend, tstart, zeroT

using IntervalArithmetic: sup
import IntervalArithmetic as IA
using IntervalBoxes: IntervalBox
using LazySets: Hyperrectangle, LazySet, overapproximate
using MathematicalSystems: AbstractContinuousSystem, BlackBoxContinuousSystem,
                           IVP, initial_state, statedim, system
using Parameters: @unpack
using ReachabilityBase.Require: @required
using Requires: @require

import ..ReachabilityAnalysis: post, numtype, rsetrep

include("FLOWSTAR.jl")
include("post.jl")

include("init.jl")

export FLOWSTAR

end  # module
