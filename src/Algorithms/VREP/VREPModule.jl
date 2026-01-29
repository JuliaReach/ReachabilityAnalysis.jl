module VREPModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, ReachSet,
                              TimeInterval, VPOLY, hasinput, zeroI, _reconvert
using ..ForwardModule: Forward

using IntervalArithmetic: (..)
using LazySets: Universe, VPolygon, VPolytope, linear_map, set
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state,
                           state_matrix, stateset
using Parameters: @unpack
using StaticArrays: SVector

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("VREP.jl")
include("post.jl")
include("reach_homog.jl")

export VREP

end  # module
