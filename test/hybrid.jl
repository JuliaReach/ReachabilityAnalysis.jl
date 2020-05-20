using ReachabilityAnalysis: _distribute, HRepIntersection,
                            WaitingList, StateInLocation, state, location

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
    prob_distr = _distribute(prob, intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr =_distribute(prob, intersection_method, check_invariant=true)
    @test !isempty(initial_state(prob_distr))

    # X0 is completely outside the invariant
    prob, _ = bouncing_ball(X0 = Hyperrectangle(low=[-10.0, -2.0], high=[-9., -1.]))
    prob_distr = _distribute(prob, intersection_method, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr = _distribute(prob, intersection_method, check_invariant=true)
    # check that this set is not added to the waiting list
    @test isempty(initial_state(prob_distr))
end

@testset "Bouncing ball: linear solvers" begin
    prob, _ = bouncing_ball()
end

@testset "Bouncing ball: nonlinear solvers" begin

end
