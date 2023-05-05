using ReachabilityAnalysis: reach_homog_QINT

@testset "QINT algorithm" begin
    # x' = x - x^2, with solution x(t) = x0*e^t / (1 - (1 - e^t)x0)
    sol = reach_homog_QINT(; a=-1.0, b=1.0, c=0.0, X0=Interval(0.1, 0.3),
                           T=10.0, Δ=0.1, δ=0.005, θ=0.01, maxiter=10_000)
    @test ρ([1.0], sol) < 1.02
end
