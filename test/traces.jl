using DifferentialEquations

@testset "DifferentialEquations solve API using an initial-value problem"
    prob = @ivp(x' = 1.0x, x(0) ∈ [1/2])

    sol = DifferentialEquations.solve(prob, tspan=(0.0, 1.0))
    @test sol isa OrdinaryDiffEq.ODECompositeSolution
    @test first(sol.t) == 0.0 && last(sol.t) == 1.0

    sol = DifferentialEquations.solve(prob, T=1.0)
    @test sol isa OrdinaryDiffEq.ODECompositeSolution
    @test first(sol.t) == 0.0 && last(sol.t) == 1.0

    # inplace (default)
    sol = DifferentialEquations.solve(prob, T=1.0, inplace=true)

    # out-of-place
    sol = DifferentialEquations.solve(prob, T=1.0, inplace=false)
end
