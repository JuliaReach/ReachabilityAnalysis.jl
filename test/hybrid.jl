using ReachabilityAnalysis: _distribute, HRepIntersection

@testset "Hybrid utility functions" begin
    intersection_method = HRepIntersection()
    # X0 belongs to the invariant
    prob, _ = bouncing_ball()
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
