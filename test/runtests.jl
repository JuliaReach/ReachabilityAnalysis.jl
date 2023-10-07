using Test
using ReachabilityAnalysis

TEST_MODELS = ["models/harmonic_oscillator.jl",
               "models/harmonic_oscillator.jl",
               "models/forced_oscillator.jl",
               "models/exponential1D.jl",
               "models/motor.jl",
               "models/linear5D.jl",
               "models/generated/VanDerPol.jl",
               "models/EMBrake.jl",
               "models/bouncing_ball.jl",
               "models/burgers.jl",
               "models/generated/Brusselator.jl",
               "models/hybrid/thermostat.jl"]

TEST_FILES = ["continuous/solve.jl",
              "continuous/symbolics.jl",
              "continuous/traces.jl",
              "discretization/discretization.jl",
              "flowpipes/flowpipes.jl",
              "flowpipes/setops.jl",
              "reachsets/reachsets.jl",
              "hybrid/hybrid.jl"]

TEST_ALGORITHMS = ["algorithms/INT.jl",
                   "algorithms/BOX.jl",
                   "algorithms/CARLIN.jl",
                   "algorithms/GLGM06.jl",
                   "algorithms/LGG09.jl",
                   "algorithms/ASB07.jl",
                   "algorithms/BFFPSV18.jl",
                   "algorithms/TMJets.jl",
                   "algorithms/ORBIT.jl",
                   "algorithms/QINT.jl",
                   "algorithms/VREP.jl",
                   "algorithms/FLOWSTAR.jl"]

foreach(include, TEST_MODELS)
foreach(include, TEST_FILES)
foreach(include, TEST_ALGORITHMS)
