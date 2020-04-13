@testset "Flowpipe interface" begin
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    sol = solve(p, tspan=(0.0, 1.0), GLGM06(δ=0.1))
    F = flowpipe(sol)

    # iteration
    @test set(F, 1) == set(F[1])
    @test set(F, length(F)) == set(F[end])

    # time span methods
    @test _isapprox(tspan(sol), 0.0 .. 1.0)
    @test tstart(F) == 0.0
    @test tend(F) ≈ 1.0
    @test tstart(F[1]) == 0.0
    @test tend(F[1]) ≈ 0.1
    @test tstart(F[end]) ≈ 0.9
    @test tend(F[end]) ≈ 1.0

    # set representation
    #@test setrep(F) isa Zonotope
    @test setrep(F) == Zonotope{Float64, Vector{Float64}, Matrix{Float64}}

    # callable behavior
    @test F(0..1) == F[1:end]
    @test F(0.05) == F[1]
    @test F(1.0) == F[end]
    @test F(0.05 .. 0.15) == F[1:2]

    # time interval not contained in flowpipe
    @test_throws ArgumentError F(0.05 .. 1.05)

    # test that in the border of the time transition we get two reachsets
    F = flowpipe(sol)
    @test F(0.1) == F[1:2]
end

@testset "Solution interface: initial states" begin
    # interval initial condition
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    solve(p, T=1.0)

    # interval initial condition
    p = @ivp(x' = -x, x(0) ∈ Interval(0, 1))
    solve(p, T=1.0)

    # deterministic initial condition, scalar for one-dimensional matrix (TODO)
    p = InitialValueProblem(@system(x' = -x), 0.5)
    solve(p, T=1.0)

    # deterministic initial condition, vector
    p = InitialValueProblem(@system(x' = -x), [0.5])
    solve(p, T=1.0)
end

@testset "Solution interface: time span" begin
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    Δt = 0.0 .. 2.0

    # only time horizon given
    sol = solve(p, T=2.0)
    @test _isapprox(tspan(sol), Δt)

    # both time horizon and time span given (TODO)
    #@test_throws ArgumentError solve(p, T=1.0, timespan=(0.0, 2.0))

    # time span given as a tuple
    sol = solve(p, tspan=(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as an interval
    sol = solve(p, tspan=0.0 .. 2.0)
    @test _isapprox(tspan(sol), Δt)

    # time span given as a LazySets interval
    sol = solve(p, tspan=Interval(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as a vector
    sol = solve(p, tspan=[0.0, 2.0])
    @test _isapprox(tspan(sol), Δt)
end

@testset "Concrete projection" begin
    prob, tspan = motor_homog()
    sol = solve(prob, tspan=tspan, GLGM06(δ=0.01), static=true)

    project(sol, (1, 3))
    project(sol, [1, 3])
    project(sol, vars=[1, 3])
    project(sol, vars=(1, 3))

    F = flowpipe(sol)
    project(F, (1, 3))
    project(F, [1, 3])
    project(F, vars=[1, 3])
    project(F, vars=(1, 3))

    R = sol[end]
    project(R, (1, 3))
    project(R, [1, 3])
    project(R, vars=(1, 3))
    project(R, vars=[1, 3])
end

# TODO:
# eachindex(F)


#=
TESTS

using LazySets, Revise, ReachabilityAnalysis
N = Float64
X = Hyperrectangle(N[0, 0], N[1, 1])
R = ReachSet(X, 0 .. 1)
set(R) ==
tstart(R) == 0.0
tend(R) == 1.0
tspan(R) == 1.0
dim(R) == 2
overapproximate(project(X, [1]), Interval) == Interval(-1, 1)
=#

@testset "Time-triggered hybrid automaton" begin
    prob = embrake_no_pv()
    s = prob.s

    # HACLD1 constructor with default zero jitter
    q = HACLD1(s.sys, s.rmap, s.Tsample)
    @test ReachabilityAnalysis.jitter(q) == 0.0

    sol = solve(prob, alg=GLGM06(δ=1e-7), T=1e-3)
    @test dim(sol) == 4
    @test flowpipe(sol) isa HybridFlowpipe
end
