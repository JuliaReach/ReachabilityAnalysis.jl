@testset "FLOWSTAR algorithm" begin
    model = "models/LotkaVolterra.model"
    ivp = @ivp(BlackBoxContinuousSystem(model, 2), x(0) ∈ (4.8 .. 5.2) × (1.8 .. 2.2))
    sol = solve(ivp, tspan=(0, 1), FLOWSTAR())
    @test sol isa ReachabilityAnalysis.ReachSolution{Flowpipe{Float64, TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}, Vector{TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}}}, FLOWSTAR}

    #=
    TODO Add example with model file generation:
    function f!(dx, x, p, t)
        dx[1] = 1.5*x[1] - x[1]*x[2]
        dx[2] = -3*x[2] + x[1]*x[2]
    end
    ivp = @ivp(x' = f!(x), dim=2, x(0) ∈ (4.8 .. 5.2) × (1.8 .. 2.2))
    sol = solve(ivp, tspan=(0, 1), δ=0.02, FLOWSTAR())
    sol isa ReachabilityAnalysis.ReachSolution{Flowpipe{Float64, TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}, Vector{TaylorModelReachSet{Float64, IntervalArithmetic.Interval{Float64}}}}, FLOWSTAR}
    =#
end