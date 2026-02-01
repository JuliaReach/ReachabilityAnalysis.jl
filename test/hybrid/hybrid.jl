using ReachabilityAnalysis: _distribute, WaitingList, StateInLocation, state, location,
                            TimeInterval, DeterministicSwitching, NonDeterministicSwitching
import Optim, OrdinaryDiffEq
using ReachabilityAnalysis: _isapprox
using LinearAlgebra: I

@testset "Hybrid utility functions" begin
    prob, _ = bouncing_ball()
    X0 = initial_state(prob)
    N = eltype(X0)
    TN = TimeInterval # TimeInterval{N}
    ST = typeof(X0)

    # converting a vector-of-tuples of the form [(X, 1), (Y, 2)] to a waiting list
    wl = convert(WaitingList{TN,ST,Int,StateInLocation{ST,Int}}, [(1, X0)])
    @test length(wl) == 1 && state(first(wl)) == X0 && location(first(wl)) == 1

    # can also pass a mixture of set types
    Z0 = convert(Zonotope, X0)
    ST = AbstractPolytope{N}
    wl = convert(WaitingList{TN,ST,Int,StateInLocation{ST,Int}}, [(1, X0), (2, Z0)])
    @test length(wl) == 2 && state(first(wl)) == X0 && location(first(wl)) == 1
    @test length(wl) == 2 && state(last(wl)) == Z0 && location(last(wl)) == 2

    # can also use [(1, X), (2, Y)]
    Z0 = convert(Zonotope, X0)
    ST = AbstractPolytope{N}
    wl = convert(WaitingList{TN,ST,Int,StateInLocation{ST,Int}}, [(X0, 1), (Z0, 2)])
    @test length(wl) == 2 && state(first(wl)) == X0 && location(first(wl)) == 1
    @test length(wl) == 2 && state(last(wl)) == Z0 && location(last(wl)) == 2

    intersection_method = HRepIntersection()

    # X0 belongs to the invariant
    prob_distr = _distribute(prob; intersection_method=intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr = _distribute(prob; intersection_method=intersection_method, check_invariant=true)
    @test !isempty(initial_state(prob_distr))

    # X0 is completely outside the invariant
    prob, _ = bouncing_ball(; X0=Hyperrectangle(; low=[-10.0, -2.0], high=[-9.0, -1.0]))
    prob_distr = _distribute(prob; intersection_method=intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr = _distribute(prob; intersection_method=intersection_method, check_invariant=true)
    # check that this set is not added to the waiting list
    @test isempty(initial_state(prob_distr))
end

@testset "Bouncing ball: linear solvers" begin
    prob, _ = bouncing_ball()
    X0 = initial_state(prob)

    # default options
    sol = solve(prob; T=3.0)
    @test setrep(sol) == HPolytope{Float64,Vector{Float64}}

    # default algorithm, but without intersecting the flowpipe with source invariant
    sol = solve(prob; T=3.0, intersect_source_invariant=false)
    @test setrep(sol) == Zonotope{Float64,Vector{Float64},Matrix{Float64}}

    # give a distribution of the initial states among the locations
    prob_a = IVP(system(prob), [(1, X0)])
    prob_b = IVP(system(prob), [(X0, 1)]) # equivalent
    sol_a = solve(prob_a; T=3.0)
    sol_b = solve(prob_b; T=3.0)

    # get vector of locations of each component flowpipe
    @test location.(sol) == [1, 1]

    # using GLGM06 + template hull intersection
    dirs = PolarDirections(10)
    sol = solve(prob; T=3.0,
                alg=GLGM06(; δ=1e-2, max_order=10, approx_model=Forward()),
                intersection_method=TemplateHullIntersection(dirs),
                clustering_method=ZonotopeClustering(),
                intersect_source_invariant=true)
end

@testset "Bouncing ball: nonlinear solvers" begin
    prob, _ = bouncing_ball()

    sol = solve(prob; T=3.0, max_jumps=1, alg=TMJets(),
                intersect_source_invariant=false,
                disjointness_method=ZonotopeEnclosure())
    @test rsetrep(sol) <: TaylorModelReachSet

    sol = solve(prob; T=3.0, max_jumps=1, alg=TMJets(),
                intersect_source_invariant=true,
                disjointness_method=ZonotopeEnclosure())
    @test setrep(sol) <: HPolytope
end

@testset "Inclusion checks" begin
    prob, dt = bouncing_ball()
    sol = solve(prob; T=5.0, alg=BOX(; δ=1e-3))

    @test sol[1][1] ⊆ HalfSpace([1.0, 0.0], 11.0)
    @test !(sol[1][1] ⊆ HalfSpace([1.0, 0.0], 10.0))

    # flowpipe
    @test sol[1] ⊆ HalfSpace([1.0, 0.0], 11.0)
    @test !(sol[1] ⊆ HalfSpace([1.0, 0.0], 10.0))

    # hybrid flowpipe
    @test sol.F ⊆ HalfSpace([1.0, 0.0], 11.0)
    @test !(sol.F ⊆ HalfSpace([1.0, 0.0], 10.0))

    # solution
    @test sol ⊆ HalfSpace([1.0, 0.0], 11.0)
    @test !(sol ⊆ HalfSpace([1.0, 0.0], 10.0))
end

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
    sol = solve(prob; alg=GLGM06(; δ=1e-7, approx_model=Forward()), max_jumps=2)
    @test dim(sol) == 4
    @test sol.alg.static == Val(false)
    @test flowpipe(sol) isa HybridFlowpipe

    sol = solve(prob;
                alg=GLGM06(; δ=1e-7, approx_model=Forward(), static=true,
                           max_order=1, dim=4, ngens=4), max_jumps=2)
    @test sol.alg.static == Val(true)

    # scenario with parameter variation
    # tested in test/algorithms/ASB07.jl
end

@ts @testset "Thermostat simulations" begin
    S = thermostat()
    T = 5.0

    # solve with default options (no simulations)
    sol = solve(S; T=T)

    # handle invariants discretely
    solsd = solve(S; T=T, ensemble=true, trajectories=10, use_discrete_callback=true)
    @test length(ensemble(solsd)) == 10

    # include the X0 vertices (here: singleton X0, so all simulations are the same)
    solsd = solve(S; T=T, ensemble=true, trajectories=10, include_vertices=true,
                  use_discrete_callback=true)
    @test length(ensemble(solsd)) == 10 + 1

    # handle invariants continuously
    solscstep = solve(S; T=T, ensemble=true, trajectories=10, dtmax=0.1)
    @test length(ensemble(solscstep)) == 10
end
