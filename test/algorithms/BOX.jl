@testset "BOX algorithm" begin

    # one-dimensional
    prob, tspan = exponential_1d()
    sol = solve(prob, tspan=tspan, BOX(δ=0.01))
    @test isa(sol.alg, BOX)
    @test setrep(sol) <: Hyperrectangle
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 1

    # higher-dimensional homogeneous
    prob, tspan = linear5D_homog()
    sol = solve(prob, tspan=tspan, BOX(δ=0.01))
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 5

    # static option
    sol = solve(prob, tspan=tspan, BOX(δ=0.01), static=true)
    HS = Hyperrectangle{Float64,StaticArrays.SArray{Tuple{5},Float64,1,5},StaticArrays.SArray{Tuple{5},Float64,1,5}}
    @test setrep(sol) == HS

    # TODO higher-dimensional with input
    #prob, tspan = linear5D()
    #sol = solve(prob, tspan=tspan, BOX(δ=0.01))
end
