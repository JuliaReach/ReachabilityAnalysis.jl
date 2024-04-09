import OrdinaryDiffEq
using ReachabilityAnalysis: _solve_ensemble, ReachSolution

@testset "Simulation traces using OrdinaryDiffEq.jl solvers" begin
    prob = @ivp(x' = 1.0x, x(0) ∈ [1 / 2])

    # `solve` extends CommonSolve.jl
    sol = OrdinaryDiffEq.solve(prob; tspan=(0.0, 1.0))
    @test sol isa ReachSolution

    # ensemble problems API
    sol = _solve_ensemble(prob; tspan=(0.0, 1.0))
    @test sol isa OrdinaryDiffEq.EnsembleSolution
    sol = _solve_ensemble(prob; T=1.0)
    @test sol isa OrdinaryDiffEq.EnsembleSolution

    # inplace (default)
    sol = _solve_ensemble(prob; T=1.0, inplace=true)

    # out-of-place
    sol = _solve_ensemble(prob; T=1.0, inplace=false)
end
