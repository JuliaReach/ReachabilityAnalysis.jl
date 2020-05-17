@testset "Hybrid utility functions" begin
    # X0 belongs to the invariant
    prob, _ = bouncing_ball()
    prob_distr = RA._distribute(prob, check_invariant=false);
    @test !isempty(initial_state(prob_distr))
    prob_distr = RA._distribute(prob, check_invariant=true)
    @test !isempty(initial_state(prob_distr))

    # X0 is completely outside the invariant
    prob, _ = bouncing_ball(X0 = Hyperrectangle(low=[-10.0, -2.0], high=[-9., -1.]))
    prob_distr = RA._distribute(prob, check_invariant=false)
    @test !isempty(initial_state(prob_distr))
    prob_distr = RA._distribute(prob, check_invariant=true)
    # check that this set is not added to the waiting list
    @test isempty(initial_state(prob_distr))
end
