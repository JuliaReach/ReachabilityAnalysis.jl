@testset "GLGM06 algorithm: homogeneous" begin

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

    # static option
    sol = solve(prob, T=5.0, GLGM06(δ=0.01, static=true))
    ZS = Zonotope{Float64, SArray{Tuple{5},Float64,1,5}, SArray{Tuple{5,5},Float64,2,25}}
    @test setrep(sol) == ZS

    # use approx model for "discrete-time" reachability
    sol = solve(prob, T=5.0, GLGM06(δ=0.01, approx_model=NoBloating()))

end

@testset "GLGM06 algorithm: inhomogeneous" begin
    # five-dim inhomogeneous
    prob, _ = linear5D()
    sol = solve(prob, T=5.0, GLGM06(δ=0.01))
    @test dim(sol) == 5
    @test tspan(shift(sol, 1.0)) == tspan(sol) + 1.0

    # use approx model for "discrete-time" reachability
    sol = solve(prob, T=5.0, GLGM06(δ=0.01, approx_model=NoBloating()))

    # motor model (8 dim)
    prob, _ = motor()

    # dense array case
    alg = GLGM06(δ=1e-2)
    sol = solve(prob, T=20.0, alg=alg)
    RT = ReachSet{Float64, Zonotope{Float64, Vector{Float64}, Matrix{Float64}}}
    @test rsetrep(sol) == RT

    # static case
    alg = GLGM06(δ=1e-2, static=true, dim=8, max_order=1, ngens=8)
    sol = solve(prob, T=20.0, alg=alg)
    RT = ReachSet{Float64, Zonotope{Float64, SVector{8, Float64}, SMatrix{8, 8, Float64, 64}}}
    @test rsetrep(sol) == RT
end
