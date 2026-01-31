using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
using StaticArrays: SVector

@testset "BOX struct" begin
    # full constructor
    alg = BOX(0.01, nothing, nothing, nothing, nothing)

    # constructor with default values
    alg = BOX(; δ=0.01)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    @test rsetrep(alg) == ReachSet{Float64,setrep(alg)}

    # StaticArrays
    alg = BOX(; δ=0.01, static=true, dim=1)
    @test setrep(alg) == Hyperrectangle{Float64,SVector{1,Float64},SVector{1,Float64}}
    @test rsetrep(alg) == ReachSet{Float64,setrep(alg)}

    # invalid construction: `static=true` requires `dim`
    # TODO should this be detected by the constructor?
    alg = BOX(; δ=0.01, static=true)
    @test_throws MethodError setrep(alg)
end

@testset "BOX algorithm" begin
    # one-dimensional
    ivp, _ = exponential_1d()
    δ = 0.01
    alg = BOX(; δ=δ)
    T = 2 * δ

    # continuous algorithm
    sol = solve(ivp; T=T, alg=alg)
    @test isa(sol.alg, BOX)
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 1
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # inhomogeneous system
    A = hcat(state_matrix(ivp))
    X = stateset(ivp)
    X0 = initial_state(ivp)
    B = Id(1, 1.0)
    U = ConstantInput(Singleton([1.0]))
    ivp2 = @ivp(ConstrainedLinearControlContinuousSystem(A, B, X, U), X(0) ∈ X0)
    sol = solve(ivp2; T=T, alg=alg)
    @test isa(sol.alg, BOX)
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 1

    # homogenization
    sol = solve(ivp2; T=T, alg=alg, homogenize=true)
    @test isa(sol.alg, BOX)
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 2
    # multiple of identity
    B = Id(1, 2.0)
    ivp3 = @ivp(ConstrainedLinearControlContinuousSystem(A, B, X, U), X(0) ∈ X0)
    sol = solve(ivp3; T=T, alg=alg, homogenize=true)
    @test isa(sol.alg, BOX)
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 2

    # invariant
    X = BallInf(zeros(1), 100.0)
    ivp3 = @ivp(ConstrainedLinearContinuousSystem(A, X), X(0) ∈ X0)
    sol = solve(ivp3; T=T, alg=alg)
    ivp3 = @ivp(ConstrainedLinearControlContinuousSystem(A, B, X, U), X(0) ∈ X0)
    sol = solve(ivp3; T=T, alg=alg)

    # recursive option
    alg = BOX(; δ=δ, recursive=true)
    sol = solve(ivp; T=T, alg=alg)
    @test_broken solve(ivp2; T=T, alg=alg)  # TODO support recursive option for inhomogeneous case

    # higher-dimensional homogeneous
    ivp, _ = linear5D_homog()
    sol = solve(ivp; T=T, alg=BOX(; δ=δ))
    @test setrep(sol) == Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}
    @test dim(sol) == 5

    # static option
    sol = solve(ivp; T=T, alg=BOX(; δ=δ), static=false) # TODO add static = true
    #HS = Hyperrectangle{Float64,StaticArrays.SArray{Tuple{5},Float64,1,5},StaticArrays.SArray{Tuple{5},Float64,1,5}} # TODO add
    #@test setrep(sol) == HS # TODO add

    # TODO higher-dimensional with input
    #ivp, _ = linear5D()
    #sol = solve(ivp, T=T, alg=BOX(δ=δ))
end
