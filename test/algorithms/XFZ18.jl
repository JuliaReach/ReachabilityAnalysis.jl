using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA

@testset "XFZ18 struct" begin
    # full constructor
    alg = RA.XFZ18Module.XFZ18{Float64,Nothing}()  # TODO remove module access once integrated

    # constructor with default values
    @test_throws MethodError RA.XFZ18Module.XFZ18()

    # struct getters
    @test RA.numtype(alg) == Float64
end

@testset "XFZ18 algorithm" begin
    alg = RA.XFZ18Module.XFZ18{Float64,Nothing}()

    # homogeneous problem
    ivp, tspan = exponential_1d()
    @test_throws MethodError solve(ivp; tspan=tspan, alg=alg)

    # inhomogeneous problem
    ivp, tspan = projectile()
    @test_throws MethodError solve(ivp; tspan=tspan, alg=alg)
end
