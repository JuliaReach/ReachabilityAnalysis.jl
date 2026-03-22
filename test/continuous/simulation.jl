using ReachabilityAnalysis, Test
using ReachabilityAnalysis: RungeKutta

@testset "Simulation" begin
    ivp, _ = exponential_1d()
    Δt = (0.0, 0.01)
    sol = solve(ivp; tspan=Δt, alg=BOX(; δ=0.01), ensemble=true)
    sim = sol.ext[:ensemble];
    @test length(sim) == 10
    sol = solve(ivp; tspan=Δt, alg=BOX(; δ=0.01), ensemble=true,
                trajectories_alg=RungeKutta(), trajectories=3)
    sim2 = sol.ext[:ensemble];
    @test length(sim) == 3
end
