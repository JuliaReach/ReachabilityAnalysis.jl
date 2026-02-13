module INTModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, ReachSet,
                              TimeInterval, hasinput, zeroT, _isdisjoint
using ..ForwardModule: Forward

using IntervalArithmetic.Symbols: (..)
using LazySets: Interval, LazySet, Universe, overapproximate, set
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, statedim, stateset
using Parameters: @unpack

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("INT.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export INT

end  # module
