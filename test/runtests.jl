using Test, ReachabilityAnalysis

# auxiliary code to skip expensive tests
begin
    __test_short = haskey(ENV, "JULIA_PKGEVAL")

   "Use the macro `ts` for tests which are deemed too slow and make the pipeline fail when running PkgEval."
    macro ts(arg)
        if !__test_short
            quote
                $(esc(arg))
            end
        end
    end

    fn(vec) = [e for e in vec if !isnothing(e)]
end

TEST_MODELS = fn(["models/harmonic_oscillator.jl",
                  "models/forced_oscillator.jl",
                  "models/exponential1D.jl",
                  "models/motor.jl",
                  "models/linear5D.jl",
                  "models/generated/VanDerPol.jl",
                  "models/EMBrake.jl",
                  "models/bouncing_ball.jl",
                  @ts("models/burgers.jl"),
                  @ts("models/generated/Brusselator.jl"),
                  @ts("models/hybrid/thermostat.jl")])

TEST_FILES = fn(["continuous/solve.jl",
                 @ts("continuous/symbolics.jl"),
                 @ts("continuous/traces.jl"),
                 @ts("discretization/exponentiation.jl"),
                 "flowpipes/flowpipes.jl",
                 @ts("flowpipes/setops.jl"),
                 "reachsets/reachsets.jl",
                 "hybrid/hybrid.jl"])

TEST_ALGORITHMS = fn(["algorithms/INT.jl",
                      "algorithms/BOX.jl",
                      @ts("algorithms/CARLIN.jl"),
                      "algorithms/GLGM06.jl",
                      "algorithms/LGG09.jl",
                      @ts("algorithms/ASB07.jl"),
                      "algorithms/BFFPSV18.jl",
                      "algorithms/TMJets.jl",
                      "algorithms/ORBIT.jl",
                      "algorithms/QINT.jl",
                      @ts("algorithms/VREP.jl"),
                      @ts("algorithms/FLOWSTAR.jl")])

foreach(include, TEST_MODELS)
foreach(include, TEST_FILES)
foreach(include, TEST_ALGORITHMS)

@ts include("Aqua.jl")
