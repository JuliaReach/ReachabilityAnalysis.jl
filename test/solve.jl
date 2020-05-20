@testset "Flowpipe interface" begin
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    δ = 0.1
    sol = solve(p, tspan=(0.0, 1.0), GLGM06(δ=δ))
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

    # time span for flowpipe slices
    @test tstart(F[1:3]) ≈ 0.0
    @test tend(F[1:3]) ≈ 3δ
    @test tspan(F[1:3]) ≈ 0.0 .. 3δ

    @test tstart(F, 1:3) ≈ 0.0  # same, more efficient
    @test tend(F, 1:3) ≈ 3δ
    @test tspan(F, 1:3) ≈ 0.0 .. 3δ

    @test tstart(sol[1:3]) ≈ 0.0
    @test tend(sol[1:3]) ≈ 3δ
    @test tspan(sol[1:3]) ≈ 0.0 .. 3δ

    # set representation
    #@test setrep(F) isa Zonotope
    N = Float64
    ZT = Zonotope{N, Vector{N}, Matrix{N}}
    @test setrep(F) == ZT
    @test rsetrep(F) == ReachSet{N, ZT}

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

    # dimension
    @test dim(F) == 1

    # different ways to input the time span or the time horizon
    sol = solve(p, tspan=(0.0, 1.0))
    sol = solve(p, (0.0, 1.0))
    sol = solve(p, tspan=[0.0, 1.0])
    sol = solve(p, [0.0, 1.0])
    sol = solve(p, T=1.0)
    sol = solve(p, 1.0)
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
    sol = solve(prob, tspan=tspan, GLGM06(δ=0.01), static=false) # TODO: add static=true as well

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

@testset "Time-triggered hybrid automaton (HACLD1)" begin
    # HACLD1 constructors
    idsys = @system(x' = Matrix(1.0I, 2, 2) * x)
    idmap = x -> x
    Ts = 1e-3

    # jitter is not defined => no jitter
    A = HACLD1(idsys, idmap, Ts)
    @inferred HACLD1(idsys, idmap, Ts)
    @test jitter(A) == 0.0 .. 0.0
    @test switching(A) == DeterministicSwitching

    # jitter is zero => no jitter
    A = HACLD1(idsys, idmap, Ts, 0.0)
    @test_broken @inferred HACLD1(idsys, idmap, Ts, 0.0)
    @test jitter(A) == 0.0 .. 0.0
    @test switching(A) == DeterministicSwitching

    # jitter is a number => symmetric jitter, [-ζ, ζ]
    A = HACLD1(idsys, idmap, Ts, 1e-8)
    @test_broken @inferred HACLD1(idsys, idmap, Ts, 1e-8)
    @test jitter(A) == -1e-8 .. 1e-8
    @test switching(A) == NonDeterministicSwitching

    # jitter is an interval [ζ⁻, ζ⁺]
    A = HACLD1(idsys, idmap, Ts, -1e-8 .. 1e-7)
    @inferred HACLD1(idsys, idmap, Ts, -1e-8 .. 1e-7)
    @test jitter(A) == -1e-8 .. 1e-7
    @test switching(A) == NonDeterministicSwitching

    # jitter is a vector [ζ⁻, ζ⁺], it is converted to an interval
    A = HACLD1(idsys, idmap, Ts, [-1e-8, 1e-7])
    @test_broken @inferred HACLD1(idsys, idmap, Ts, [-1e-8, 1e-7])
    @test _isapprox(jitter(A), -1e-8 .. 1e-7)
    @test switching(A) == NonDeterministicSwitching

    # jitter is a tuple (ζ⁻, ζ⁺), it is converted to an interval
    A = HACLD1(idsys, idmap, Ts, (-1e-8, 1e-7))
    @test_broken @inferred HACLD1(idsys, idmap, Ts, (-1e-8, 1e-7))
    @test _isapprox(jitter(A), -1e-8 .. 1e-7)
    @test switching(A) == NonDeterministicSwitching
end

@testset "Time-triggered solve (EMBrake)" begin
    # scenario without parameter variation
    prob = embrake_no_pv()
    sol = solve(prob, alg=GLGM06(δ=1e-7), max_jumps=1)
    @test dim(sol) == 4
    @test flowpipe(sol) isa HybridFlowpipe

    # scenario with parameter variation
    # tested in test/algorithms/ASB07.jl
end
