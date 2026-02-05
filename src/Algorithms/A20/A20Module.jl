module A20Module

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, NoEnclosure,
                              ReachSet, TimeInterval, hasinput, zeroT,
                              _reconvert
using ..FirstOrderZonotopeModule: FirstOrderZonotope
using ..Overapproximate: _convert_or_overapproximate
using LazySets: GIR05, LazySet, Zonotope, reduce_order
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, stateset
using Parameters: @unpack

import ..ReachabilityAnalysis: post, numtype, rsetrep, setrep, step_size

include("A20.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export A20

end  # module
