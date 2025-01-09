@testset "INT algorithm" begin
    ivp, tspan = exponential_1d()
    alg = INT(; δ=0.01)
    # continuous algorithm
    sol = solve(ivp; tspan=tspan, alg=alg)
    @test isa(sol.alg, INT)
    @test dim(sol) == 1
    @test length(sol) == 100
    @test setrep(sol) <: Interval
    @test setrep(sol) == Interval{Float64}
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    ivp, tspan = exponential_1d(; invariant=HalfSpace([-1.0], -0.3)) # x >= 0.3
    sol_inv = solve(ivp; tspan=tspan, alg=INT(; δ=0.01))
    @test length(sol_inv) == 52
    @test [0.3] ∈ sol_inv[end]
    # check that the following reach-set escapes the invariant
    @test [0.3] ∈ sol[52] && [0.3] ∉ sol[53]

    # check NSTEPS option
    sol = solve(ivp; NSTEPS=10, alg=INT(; δ=0.01))
    @test length(sol) == 10

    # doesn't work for higher dimensional systems
    ivp, tspan = linear5D_homog()
    @test_throws ArgumentError solve(ivp, tspan=tspan, alg=INT(; δ=0.01))
end

# Issue #378
@testset "Scalar affine ODE" begin
    ivp = @ivp(x' = x + [1.0], x(0) ∈ Singleton([1.0]))
    sol = solve(ivp; T=4.0)
end
