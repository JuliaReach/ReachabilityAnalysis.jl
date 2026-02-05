using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA

@testset "QINT struct" begin
    # full constructor
    alg = QINT(0.0, 0.01, 0.0, 0, nothing)

    # constructor with default values
    alg = QINT(; Δ=1.0, δ=0.01, θ=1.0, maxiter=10)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test setrep(alg) == Interval{Float64}
    @test rsetrep(alg) == ReachSet{Float64,Interval{Float64}}
end

@testset "QINT algorithm" begin
    ivp, tspan = exponential_1d()
    alg = QINT(; Δ=0.1, δ=0.01, θ=0.01, maxiter=10)
    @test_broken RA.reach_homog_QINT(alg, ivp, tspan) isa HybridFlowpipe  # TODO finish implementation; method should also be called `post` so one can call `solve`

    # call `reach_homog_QINT` directly
    # x' = x - x^2, with solution x(t) = x0*e^t / (1 - (1 - e^t)x0)
    sol = RA.reach_homog_QINT(; a=-1.0, b=1.0, c=0.0, X0=Interval(0.1, 0.3),
                              T=10.0, Δ=0.1, δ=0.005, θ=0.01, maxiter=10_000)
    @test ρ([1.0], sol) < 1.02
end
