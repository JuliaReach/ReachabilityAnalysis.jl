using ReachabilityAnalysis, Test
using ReachabilityAnalysis: RungeKutta

@testset "Simulation" begin
    ivp, _ = exponential_1d()
    Δt = (0.0, 0.1)
    sol = solve(ivp; tspan=Δt, alg=BOX(δ=0.01), ensemble=true, trajectories_alg=RungeKutta(4, 1e-3))
    sim = sol.ext[:ensemble]
    @test length(sim) == 10
end
