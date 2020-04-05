@testset "INT algorithm" begin
    prob, tspan = exponential_1d()
    sol = solve(prob, tspan=tspan, INT(δ=0.01))
    @test isa(sol.alg, INT)
    @test setrep(sol) <: Interval
    @test setrep(sol) == Interval{Float64,IntervalArithmetic.Interval{Float64}}

    # doesn't work for higher dimensional systems
    prob, tspan = linear5D_homog()
    @test_throws ArgumentError solve(prob, tspan=tspan, INT(δ=0.01))
end
