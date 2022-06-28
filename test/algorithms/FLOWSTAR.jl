@testset "FLOWSTAR algorithm" begin
    # Passing a model file.
    model = joinpath(@__DIR__, "..", "models", "LotkaVolterra.model")
    ivp = @ivp(BlackBoxContinuousSystem(model, 2), x(0) ∈ (4.8 .. 5.2) × (1.8 .. 2.2))
    sol = solve(ivp, tspan=(0, 1), FLOWSTAR(δ=0.02))
    RT = TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}
    @test sol isa ReachabilityAnalysis.ReachSolution{Flowpipe{Float64, RT, Vector{RT}},
                    FLOWSTAR{Float64, Flowstar.FixedTMOrder, Flowstar.QRPreconditioner, Flowstar.NonPolyODEScheme}}

    # Generating the model file from the given system.
    function f!(dx, x, p, t)
        dx[1] = 1.5*x[1] - x[1]*x[2]
        dx[2] = -3*x[2] + x[1]*x[2]
    end
    ivp = @ivp(x' = f!(x), dim=2, x(0) ∈ (4.8 .. 5.2) × (1.8 .. 2.2))
    sol = solve(ivp, tspan=(0, 1), FLOWSTAR(δ=0.02))
    RT = TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}
    @test sol isa ReachabilityAnalysis.ReachSolution{Flowpipe{Float64, RT, Vector{RT}},
                    FLOWSTAR{Float64, Flowstar.FixedTMOrder, Flowstar.QRPreconditioner, Flowstar.NonPolyODEScheme}}
end
