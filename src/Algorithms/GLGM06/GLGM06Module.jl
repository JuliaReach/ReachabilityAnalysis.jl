module GLGM06Module

using ..ReachabilityAnalysis: AbstractContinuousPost, AbstractDisjointnessMethod,
                              Flowpipe, NoEnclosure, ReachSet, TimeInterval,
                              hasinput, linear_map_minkowski_sum, zeroT,
                              _isdisjoint, _reconvert
using ..FirstOrderZonotopeModule: FirstOrderZonotope
using ..Overapproximate: _convert_or_overapproximate

using IntervalArithmetic: (..)
using LazySets: AbstractReductionMethod, GIR05, LazySet, Universe, Zonotope,
                linear_map, linear_map!, reduce_order, set
using LinearAlgebra: mul!
using MathematicalSystems: AbstractDiscreteSystem, IVP, initial_state, inputset,
                           state_matrix, stateset
using Parameters: @unpack
using StaticArrays: SMatrix, SVector

import ..ReachabilityAnalysis: post, numtype, setrep, step_size, rsetrep

include("GLGM06.jl")
include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

export GLGM06

end  # module
