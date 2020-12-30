using Test, ReachabilityAnalysis, StaticArrays

using ReachabilityAnalysis: _isapprox, setrep, rsetrep,
                           DeterministicSwitching, NonDeterministicSwitching,
                           BoxEnclosure, ZonotopeEnclosure

using Polyhedra, CDDLib # for VREP algorithm
# fix namespace conflicts with Polyhedra
using LazySets: dim, HalfSpace, Interval, Line2D, translate, project
const RA = ReachabilityAnalysis

import IntervalArithmetic
const IA = IntervalArithmetic

# load test models
include("models/harmonic_oscillator.jl")
include("models/forced_oscillator.jl")
include("models/exponential1D.jl")
include("models/motor.jl")
include("models/linear5D.jl")
include("models/VanDerPol.jl")
include("models/EMBrake.jl")
include("models/bouncing_ball.jl")
include("models/burgers.jl")

#include("utils.jl")
include("solve.jl")
include("reachsets.jl")
include("flowpipes.jl")
#include("traces.jl")
include("setops.jl")
include("arrayops.jl")

# algorithms
include("algorithms/INT.jl")
include("algorithms/BOX.jl")
include("algorithms/GLGM06.jl")
include("algorithms/LGG09.jl")
include("algorithms/ASB07.jl")
include("algorithms/BFFPSV18.jl")
include("algorithms/TMJets.jl")
include("algorithms/ORBIT.jl")
include("algorithms/QINT.jl")
include("algorithms/VREP.jl")

include("hybrid.jl")
