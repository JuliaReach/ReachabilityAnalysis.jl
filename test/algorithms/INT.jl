@testset "INT algorithm" begin
    prob, tspan = exponential_1d()
    sol = solve(prob; tspan=tspan, INT(; δ=0.01))
    @test isa(sol.alg, INT)
    @test dim(sol) == 1
    @test length(sol) == 100
    @test setrep(sol) <: Interval
    @test setrep(sol) == Interval{Float64,IntervalArithmetic.Interval{Float64}}

    prob, tspan = exponential_1d(; invariant=HalfSpace([-1.0], -0.3)) # x >= 0.3
    sol_inv = solve(prob; tspan=tspan, INT(; δ=0.01))
    @test length(sol_inv) == 52
    @test [0.3] ∈ sol_inv[end]
    # check that the following reach-set escapes the invariant
    @test [0.3] ∈ sol[52] && [0.3] ∉ sol[53]

    # check NSTEPS option
    sol = solve(prob; NSTEPS=10, INT(; δ=0.01))
    @test length(sol) == 10

    # doesn't work for higher dimensional systems
    prob, tspan = linear5D_homog()
    @test_throws ArgumentError solve(prob, tspan=tspan, INT(; δ=0.01))
end

# Issue #378
@testset "Scalar affine ODE" begin
    prob = @ivp(x' = x + [1.0], x(0) ∈ Singleton([1.0]))
    sol = solve(prob; T=4.0)
end
