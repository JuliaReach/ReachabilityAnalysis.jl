using Test, ReachabilityAnalysis, StaticArrays

using ReachabilityAnalysis: _isapprox, setrep, rsetrep,
                           DeterministicSwitching, NonDeterministicSwitching,
                           BoxEnclosure, ZonotopeEnclosure

import IntervalArithmetic
const IA = IntervalArithmetic

# load test models
include("models/exponential1D.jl")
include("models/motor.jl")
include("models/linear5D.jl")
include("models/vanderpol.jl")
include("models/embrake.jl")
include("models/bouncing_ball.jl")

#include("utils.jl")
include("solve.jl")
include("reachsets.jl")
include("flowpipes.jl")
#include("traces.jl")
include("setops.jl")

include("algorithms/INT.jl")
include("algorithms/BOX.jl")
include("algorithms/GLGM06.jl")
include("algorithms/LGG09.jl")
include("algorithms/ASB07.jl")
include("algorithms/BFFPSV18.jl")
include("algorithms/TMJets.jl")

include("hybrid.jl")
