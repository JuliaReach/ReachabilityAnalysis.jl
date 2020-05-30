@testset "GLGM06 algorithm" begin

    # one-dimensional
    prob, _ = exponential_1d()
    sol = solve(prob, T=5.0, GLGM06(δ=0.01))
    @test isa(sol.alg, GLGM06)
    @test setrep(sol) <: Zonotope
    @test setrep(sol) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}
    @test dim(sol) == 1

    # higher-dimensional homogeneous
    prob, _ = linear5D_homog()
    sol = solve(prob, T=5.0, GLGM06(δ=0.01))
    @test dim(sol) == 5

    # static option TODO
    #sol = solve(prob, tspan=tspan, GLGM06(δ=0.01), static=true)
    #ZS = Zonotope{Float64,StaticArrays.SArray{Tuple{5},Float64,1,5},StaticArrays.SArray{Tuple{5,5},Float64,2,25}}
    #@test setrep(sol) == ZS

    # TODO higher-dimensional with input
    #prob, tspan = linear5D()
    #sol = solve(prob, tspan=tspan, GLGM06(δ=0.01))

    # use approx model for "discrete-time" reachability
    sol = solve(prob, T=5.0, GLGM06(δ=0.01, approx_model=NoBloating()))

    # higher-dimensional inhomogeneous
    prob, _ = linear5D_homog()
    sol = solve(prob, T=5.0, GLGM06(δ=0.01))
    @test dim(sol) == 5

    # use approx model for "discrete-time" reachability
    sol = solve(prob, T=5.0, GLGM06(δ=0.01, approx_model=NoBloating()))
end
