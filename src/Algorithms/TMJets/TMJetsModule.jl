module TMJetsModule

using ..ReachabilityAnalysis: AbstractContinuousPost, AbstractDisjointnessMethod,
                              Flowpipe, TimeInterval, TaylorModelReachSet,
                              VectorField, ZonotopeEnclosure, domain,
                              inplace_field!, tend, tstart, zeroT, _isdisjoint,
                              _normalize, _overapproximate_structured,
                              _overapproximate_structured_full, _shrink_wrapping

import IntervalArithmetic as IA
using IntervalBoxes: IntervalBox
using LazySets: AbstractHyperrectangle, AbstractZonotope, CartesianProduct,
                Interval, LazySet, Zonotope, box_approximation, dim, order,
                overapproximate, set
using LinearAlgebra: isdiag
using MathematicalSystems: AbstractContinuousSystem, IVP, initial_state,
                           isaffine, islinear, statedim, stateset
using Parameters: @unpack, @with_kw
using TaylorModels: TaylorN, TaylorModel1, remainder, set_variables
using TaylorModels.ValidatedInteg: _DEF_MINABSTOL, validated_integ,
                                   validated_integ2

import ..ReachabilityAnalysis: post, numtype, rsetrep

include("common.jl")
include("TMJets21a/TMJets21a.jl")
include("TMJets21b/TMJets21b.jl")

"""
    TMJets

The algorithm TMJets defaults to `TMJets21b`.
"""
const TMJets = TMJets21b

end  # module
