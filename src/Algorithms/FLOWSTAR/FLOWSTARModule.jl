module FLOWSTARModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, TimeInterval,
                              TaylorModelReachSet, domain, tend, tstart, zeroT

using IntervalArithmetic: sup
import IntervalArithmetic as IA
using MathematicalSystems: AbstractContinuousSystem, BlackBoxContinuousSystem,
                           IVP
using Parameters: @unpack
using ReachabilityBase.Require: require

import ..ReachabilityAnalysis: post, numtype, rsetrep

include("FLOWSTAR.jl")
include("post.jl")

export FLOWSTAR

end  # module
