@testset "BOX algorithm" begin

    # one-dimensional
    ivp, tspan = exponential_1d()
    alg = BOX(; δ=0.01)
    # continuous algorithm
    sol = solve(ivp; tspan=tspan, alg=alg)
    @test isa(sol.alg, BOX)
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 1
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # homogenization
    A = hcat(state_matrix(ivp))
    X = stateset(ivp)
    X0 = initial_state(ivp)
    B = Id(1, 1.0)
    U = ConstantInput(Singleton([1.0]))
    ivp2 = @ivp(ConstrainedLinearControlContinuousSystem(A, B, X, U), X(0) ∈ X0)
    sol2 = solve(ivp2; tspan=(0.0, 1.0), alg=alg, homogenize=true)
    @test isa(sol2.alg, BOX)
    @test setrep(sol2) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol2) == 2
    # multiple of identity
    B = Id(1, 2.0)
    ivp3 = @ivp(ConstrainedLinearControlContinuousSystem(A, B, X, U), X(0) ∈ X0)
    sol3 = solve(ivp3; tspan=(0.0, 1.0), alg=alg, homogenize=true)
    @test isa(sol3.alg, BOX)
    @test setrep(sol3) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol3) == 2

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
