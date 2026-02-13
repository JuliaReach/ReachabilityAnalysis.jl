module CARLINModule

using ..ReachabilityAnalysis: AbstractContinuousPost, Flowpipe, LGG09,
                              MixedFlowpipe, ReachSet, flowpipe, solve, tspan,
                              zeroI, zeroT, _get_T
using ..ForwardModule: Forward

using Base: @kwdef
using CarlemanLinearization: build_matrix, error_bound_specabs
using IntervalArithmetic: inf, mince, sup
using IntervalArithmetic.Symbols: (..)
import IntervalArithmetic as IA
using IntervalBoxes: IntervalBox
using LazySets: AbstractHyperrectangle, BallInf, CartesianProductArray,
                Hyperrectangle, Interval, box_approximation, dim, high, low,
                project, set, symmetric_interval_hull, ⊕
using LazySets.Approximations: CustomDirections
using MathematicalSystems: AbstractContinuousSystem, BlackBoxContinuousSystem,
                           IVP, LinearContinuousSystem, initial_state, system
using Parameters: @unpack
using ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityBase.Require: @required
using Requires: @require

import ..ReachabilityAnalysis: post
import CarlemanLinearization: lift_vector
import MathematicalSystems: statedim

include("CARLIN.jl")
include("kronecker.jl")
include("post.jl")
include("reach.jl")

include("init.jl")

export CARLIN,
       CanonicalQuadraticForm

end  # module
