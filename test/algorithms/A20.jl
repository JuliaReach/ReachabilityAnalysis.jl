using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA

@testset "A20 struct" begin
    # full constructor
    alg = A20(0.01, nothing, 1)

    # constructor with default values
    alg = A20(; δ=0.01)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test rsetrep(alg) == ReachSet{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
end

@testset "A20 algorithm" begin
    alg = A20(; δ=0.01)

    # homogeneous problem
    ivp, tspan = exponential_1d()
    sol = solve(ivp; tspan=tspan, alg=alg)
    @test sol.alg == alg
    @test dim(sol) == 1
    @test length(sol) == 100
    @test setrep(sol) == Zonotope{Float64,Vector{Float64},Matrix{Float64}}

    # inhomogeneous problem
    ivp, tspan = projectile()
    sol = solve(ivp; tspan=tspan, alg=alg)
    @test sol.alg == alg
    @test dim(sol) == 4
    @test length(sol) == 2000
    @test setrep(sol) == Zonotope{Float64,Vector{Float64},Matrix{Float64}}
end
