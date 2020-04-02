using Test, ReachabilityAnalysis
using ReachabilityAnalysis: _isapprox

#include("utils.jl")
include("solve.jl")
include("reachsets.jl")
#include("traces.jl")

# load test models
include("models/linear/exponential1D.jl")
include("models/nonlinear/vanderpol.jl")

# check algorithms
#include("algorithms/INT.jl")
#include("algorithms/BOX.jl")
#include("algorithms/GLGM06.jl")
#include("algorithms/BFFPSV18.jl")
include("algorithms/TMJets.jl")
