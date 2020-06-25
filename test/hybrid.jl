using ReachabilityAnalysis: _distribute, WaitingList, StateInLocation, state, location

@testset "Hybrid utility functions" begin

    prob, _ = bouncing_ball()
    X0 = initial_state(prob)
    N = eltype(X0)
    ST = typeof(X0)

    # converting a vector-of-tuples of the form [(X, 1), (Y, 2)] to a waiting list
    wl = convert(WaitingList{N, ST, Int, StateInLocation{ST, Int}}, [(1, X0)])
    @test length(wl) == 1 && state(first(wl)) == X0 && location(first(wl)) == 1

    # can also pass a mixture of set types
    Z0 = convert(Zonotope, X0)
    ST = AbstractPolytope{N}
    wl = convert(WaitingList{N, ST, Int, StateInLocation{ST, Int}}, [(1, X0), (2, Z0)])
    @test length(wl) == 2 && state(first(wl)) == X0 && location(first(wl)) == 1
    @test length(wl) == 2 && state(last(wl)) == Z0 && location(last(wl)) == 2

    # can also use [(1, X), (2, Y)]
    Z0 = convert(Zonotope, X0)
    ST = AbstractPolytope{N}
    wl = convert(WaitingList{N, ST, Int, StateInLocation{ST, Int}}, [(X0, 1), (Z0, 2)])
    @test length(wl) == 2 && state(first(wl)) == X0 && location(first(wl)) == 1
    @test length(wl) == 2 && state(last(wl)) == Z0 && location(last(wl)) == 2

    intersection_method = HRepIntersection()

    # X0 belongs to the invariant
    prob_distr = _distribute(prob, intersection_method=intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr =_distribute(prob, intersection_method=intersection_method, check_invariant=true)
    @test !isempty(initial_state(prob_distr))

    # X0 is completely outside the invariant
    prob, _ = bouncing_ball(X0 = Hyperrectangle(low=[-10.0, -2.0], high=[-9., -1.]))
    prob_distr = _distribute(prob, intersection_method=intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr = _distribute(prob, intersection_method=intersection_method, check_invariant=true)
    # check that this set is not added to the waiting list
    @test isempty(initial_state(prob_distr))
end

@testset "Bouncing ball: linear solvers" begin
    prob, _ = bouncing_ball()
    X0 = initial_state(prob)

    # default options
    sol = solve(prob, T=3.0)
    @test setrep(sol) == HPolytope{Float64, Vector{Float64}}

    # default algorithm, but without intersecting the flowpipe with source invariant
    sol = solve(prob, T=3.0, intersect_source_invariant=false)
    @test setrep(sol) == Zonotope{Float64,Vector{Float64},Matrix{Float64}}

    # give a distribution of the initial states among the locations
    prob_a = IVP(system(prob), [(1, X0)])
    prob_b = IVP(system(prob), [(X0, 1)]) # equivalent
    sol_a = solve(prob_a, T=3.0)
    sol_b = solve(prob_b, T=3.0)

    # get vector of locations of each component flowpipe
    @test location.(sol) == [1, 1]

    # using GLGM06 + template hull intersection
    dirs = PolarDirections(10)
    sol = solve(prob, T=3.0, alg=GLGM06(δ=1e-2, max_order=10),
                intersection_method=TemplateHullIntersection(dirs),
                clustering_method=ZonotopeClustering(),
                intersect_source_invariant=true)
end

@testset "Bouncing ball: nonlinear solvers" begin
    prob, _ = bouncing_ball()

    sol = solve(prob, T=3.0, max_jumps=1, alg=TMJets(),
                intersect_source_invariant=false,
                disjointness_method=ZonotopeEnclosure())
    @test rsetrep(sol) <: TaylorModelReachSet

    sol = solve(prob, T=3.0, max_jumps=1, alg=TMJets(),
                intersect_source_invariant=true,
                disjointness_method=ZonotopeEnclosure())
    @test setrep(sol) <: HPolytope
end

@testset "Inclusion checks" begin
    prob, dt = bouncing_ball()
    sol = solve(prob, T=5.0, alg=BOX(δ=1e-3))

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
