@testset "BOX algorithm" begin

    # one-dimensional
    ivp, tspan = exponential_1d()
    alg = BOX(; δ=0.01)
    # continuous algorithm
    sol = solve(ivp; tspan=tspan, alg=alg)
    @test isa(sol.alg, BOX)
    @test setrep(sol) <: Hyperrectangle
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 1
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # higher-dimensional homogeneous
    ivp, tspan = linear5D_homog()
    sol = solve(ivp; tspan=tspan, alg=BOX(; δ=0.01))
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 5

    # static option
    sol = solve(ivp; tspan=tspan, alg=BOX(; δ=0.01), static=false) # TODO add static = true
    #HS = Hyperrectangle{Float64,StaticArrays.SArray{Tuple{5},Float64,1,5},StaticArrays.SArray{Tuple{5},Float64,1,5}} # TODO add
    #@test setrep(sol) == HS # TODO add

    # TODO higher-dimensional with input
    #ivp, tspan = linear5D()
    #sol = solve(ivp, tspan=tspan, alg=BOX(δ=0.01))
end
