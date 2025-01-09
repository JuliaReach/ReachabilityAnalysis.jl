@testset "GLGM06 algorithm: homogeneous" begin

    # one-dimensional
    ivp, _ = exponential_1d()
    alg = GLGM06(; δ=0.01)
    # continuous algorithm
    sol = solve(ivp; T=5.0, alg=alg)
    @test isa(sol.alg, GLGM06)
    @test setrep(sol) <: Zonotope
    @test setrep(sol) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}
    @test dim(sol) == 1
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # higher-dimensional homogeneous
    ivp, _ = linear5D_homog()
    sol = solve(ivp; T=5.0, alg=alg)
    @test dim(sol) == 5

    # static option
    @ts begin
        sol = solve(ivp; T=5.0, alg=GLGM06(; δ=0.01, static=true))
        ZS = Zonotope{Float64,SVector{5,Float64},SMatrix{5,16,Float64,80}}
        @test setrep(sol) == ZS
    end

    # use approx model for "discrete-time" reachability
    sol = solve(ivp; T=5.0, alg=GLGM06(; δ=0.01, approx_model=NoBloating()))
end

@testset "GLGM06 algorithm: inhomogeneous" begin
    # five-dim inhomogeneous
    ivp, _ = linear5D()
    sol = solve(ivp; T=5.0, alg=GLGM06(; δ=0.01))
    @test dim(sol) == 5
    @test tspan(shift(sol, 1.0)) == tspan(sol) + 1.0

    # use approx model for "discrete-time" reachability
    sol = solve(ivp; T=5.0, alg=GLGM06(; δ=0.01, approx_model=NoBloating()))

    # motor model (8 dim)
    ivp, _ = motor()

    # dense array case
    alg = GLGM06(; δ=1e-2)
    sol = solve(ivp; T=20.0, alg=alg)
    RT = ReachSet{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}
    @test rsetrep(sol) == RT

    # static case
    @ts begin
        alg = GLGM06(; δ=1e-2, static=true, dim=8, max_order=1, ngens=8)
        sol = solve(ivp; T=20.0, alg=alg)
        RT = ReachSet{Float64,Zonotope{Float64,SVector{8,Float64},SMatrix{8,8,Float64,64}}}
        @test rsetrep(sol) == RT
    end
end
