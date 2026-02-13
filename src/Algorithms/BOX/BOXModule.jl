module BOXModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, ReachSet,
                              TimeInterval, hasinput, zeroT, _isdisjoint,
                              _reconvert
using ..ForwardModule: Forward
using ..Overapproximate: _overapproximate

using IntervalArithmetic: (..)
using LazySets: Hyperrectangle, LazySet, Universe, dim, overapproximate, set
using LinearAlgebra: mul!
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, stateset
using Parameters: @unpack
using StaticArrays: SVector

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("BOX.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export BOX

end  # module
