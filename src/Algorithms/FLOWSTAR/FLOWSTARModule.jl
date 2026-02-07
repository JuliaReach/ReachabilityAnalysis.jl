module FLOWSTARModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, TimeInterval,
                              TaylorModelReachSet, domain, tend, tstart, zeroI

using IntervalArithmetic: IntervalBox
import IntervalArithmetic as IA
using LazySets: Hyperrectangle, overapproximate
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
