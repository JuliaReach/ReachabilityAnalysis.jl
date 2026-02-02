using ReachabilityAnalysis, Test

@testset "Homogenization" begin
    # AffineContinuousSystem
    A = [0.0 1; -1 0]
    b = [0.1, 0.2]
    X0 = BallInf(zeros(2), 1.0)
    s = @system(x' = A * x + b)
    ivp = @ivp(s, x(0) ∈ X0)
    s2 = homogenize(s)
    @test s2 isa LinearContinuousSystem && state_matrix(s2) == [A b; zeros(1, 3)]
    ivp2 = homogenize(ivp)
    @test system(ivp2) == s2 && initial_state(ivp2) == X0 × Singleton([1.0])
end
