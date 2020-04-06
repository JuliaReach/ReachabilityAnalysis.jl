using Test, ReachabilityAnalysis, StaticArrays
using ReachabilityAnalysis: _isapprox

# load test models
include("models/linear/exponential1D.jl")
include("models/linear/motor.jl")
include("models/linear/linear5D.jl")
include("models/nonlinear/vanderpol.jl")

#include("utils.jl")
include("solve.jl")
include("reachsets.jl")
#include("traces.jl")

include("algorithms/INT.jl")
include("algorithms/BOX.jl")
include("algorithms/GLGM06.jl")
#include("algorithms/ASB07.jl")
#include("algorithms/BFFPSV18.jl")
include("algorithms/TMJets.jl")
