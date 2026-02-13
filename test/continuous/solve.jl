using StaticArrays: SArray
using ReachabilityAnalysis: TimeInterval, _isapprox

@testset "Default continuous post-operator" begin
    prob, _ = motor_homog()
    sol = solve(prob; T=0.001)
    @test sol.alg isa GLGM06
    @test sol.alg.static == Val{false}()
    @test setrep(sol) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}

    # static case
    sol = solve(prob; T=0.001, static=true)
    @test sol.alg isa GLGM06
    @test sol.alg.static == Val{true}()
    @test setrep(sol) ==
          Zonotope{Float64,SArray{Tuple{8},Float64,1,8},SArray{Tuple{8,13},Float64,2,104}}

    # if the static kwarg is passed outside the algorithm => it is ignored
    sol = solve(prob; T=0.001, alg=GLGM06(; δ=1e-5), static=true)
    @test setrep(sol) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}
end

@testset "Flowpipe interface" begin
    p = @ivp(x' = -x, x(0) ∈ Interval(0, 1))
    δ = 0.1
    sol = solve(p; tspan=(0.0, 1.0), alg=GLGM06(; δ=δ))
    F = flowpipe(sol)

    # iteration
    @test set(F, 1) == set(F[1])
    @test set(F, length(F)) == set(F[end])

    # time span methods
    @test _isapprox(tspan(sol), TimeInterval(0.0 .. 1.0))
    @test tstart(F) == 0.0
    @test tend(F) ≈ 1.0
    @test tstart(F[1]) == 0.0
    @test tend(F[1]) ≈ 0.1
    @test tstart(F[end]) ≈ 0.9
    @test tend(F[end]) ≈ 1.0

    # time span for flowpipe slices
    @test tstart(F[1:3]) ≈ 0.0
    @test tend(F[1:3]) ≈ 3δ
    @test _isapprox(tspan(F[1:3]), TimeInterval(0.0 .. 3δ))

    @test tstart(F, 1:3) ≈ 0.0  # same, more efficient
    @test tend(F, 1:3) ≈ 3δ
    @test _isapprox(tspan(F, 1:3), TimeInterval(0.0 .. 3δ))

    @test tstart(sol[1:3]) ≈ 0.0
    @test tend(sol[1:3]) ≈ 3δ
    @test _isapprox(tspan(sol[1:3]), TimeInterval(0.0 .. 3δ))

    # set representation
    @test setrep(F) <: Zonotope
    N = Float64
    ZT = Zonotope{N,Vector{N},Matrix{N}}
    @test setrep(F) == ZT
    @test rsetrep(F) == ReachSet{N,ZT}

    # callable behavior
    @test F(0 .. 1) == F[1:end]
    @test F(0.05) == F[1]
    @test F(1.0) == F[end]
    @test F(0.05 .. 0.15) == F[1:2]

    # time interval *not* strictly contained in flowpipe returns a subset
    @test tspan(F(0.05 .. 1.05)) == tspan(F)

    # time interval totally outside the flowpipe returns an error
    @test_throws ArgumentError F(2 .. 3)

    # test that in the border of the time transition we get two reachsets
    F = flowpipe(sol)
    @test F(0.1) == F[1:2]

    # dimension
    @test dim(F) == 1

    # different ways to input the time span or the time horizon
    sol = solve(p; tspan=(0.0, 1.0))
    sol = solve(p, (0.0, 1.0))
    sol = solve(p; tspan=[0.0, 1.0])
    sol = solve(p, [0.0, 1.0])
    sol = solve(p; T=1.0)
    sol = solve(p, 1.0)
end

@testset "Solution interface: initial states" begin
    # interval initial condition
    p = @ivp(x' = -x, x(0) ∈ Interval(0, 1))
    solve(p; T=1.0)

    # deterministic initial condition
    p = InitialValueProblem(@system(x' = -x), Singleton([0.5]))
    solve(p; T=1.0)
end

@testset "Solution interface: time span" begin
    p = @ivp(x' = -x, x(0) ∈ Interval(0, 1))
    Δt = TimeInterval(0.0 .. 2.0)

    # only time horizon given
    sol = solve(p; T=2.0)
    @test _isapprox(tspan(sol), Δt)

    # both time horizon and time span given => error
    @test_throws ArgumentError solve(p, T=1.0, tspan=(0.0, 2.0))

    # time span given as a tuple
    sol = solve(p; tspan=(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as an interval
    sol = solve(p; tspan=0.0 .. 2.0)
    @test _isapprox(tspan(sol), Δt)

    # time span given as a LazySets interval
    sol = solve(p; tspan=Interval(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as a vector
    sol = solve(p; tspan=[0.0, 2.0])
    @test _isapprox(tspan(sol), Δt)
end

@testset "Concrete projection" begin
    prob, tspan = motor_homog()
    sol = solve(prob; T=0.001, alg=GLGM06(; δ=1e-5), static=false)

    project(sol, (1, 3))
    project(sol, [1, 3])
    project(sol; vars=[1, 3])
    project(sol; vars=(1, 3))

    F = flowpipe(sol)
    project(F, (1, 3))
    project(F, [1, 3])
    project(F; vars=[1, 3])
    project(F; vars=(1, 3))

    R = sol[end]
    project(R, (1, 3))
    project(R, [1, 3])
    project(R; vars=(1, 3))
    project(R; vars=[1, 3])

    # TODO concrete projection of sets with static array fields
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
