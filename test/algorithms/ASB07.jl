using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
import ReachabilityAnalysis.Exponentiation as EXP
using StaticArrays: SMatrix, SVector

@testset "ASB07 struct" begin
    # full constructor
    alg = ASB07(0.01, nothing, 1, nothing, nothing, nothing, nothing, nothing)

    # constructor with default values
    alg = ASB07(; δ=0.01)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test rsetrep(alg) == ReachSet{Float64,Zonotope{Float64,Vector{Float64},Matrix{Float64}}}

    # StaticArrays
    alg = ASB07(; δ=0.01, static=true, dim=1, ngens=1)
    @test rsetrep(alg) ==
          ReachSet{Float64,Zonotope{Float64,SVector{1,Float64},SMatrix{1,1,Float64,1}}}

    # invalid construction: `static=true` requires `dim` and `ngens`
    # TODO should this be detected by the constructor?
    alg = ASB07(; δ=0.01, static=true)
    @test_throws MethodError rsetrep(alg)
end

@testset "ASB07 algorithm: homogeneous case" begin
    # [AlthoffSB07; Example 1](@citet)

    # linear ODE: x' = Ax
    A = IntervalMatrix([-1.0±0.05 -4.0±0.05;
                        4.0±0.05 -1.0±0.05])
    X0 = BallInf([1.0, 1.0], 0.1)
    ivp = @ivp(x' = A * x, x(0) ∈ X0)
    alg = ASB07(; δ=0.04)

    # continuous algorithm
    sol1 = solve(ivp; tspan=(0.0, 1.0), alg=alg)
    @test length(sol1) == 25
    @test sol1.alg.recursive == Val(true) # default is recursive TODO add isrecursive ?
    @test dim(sol1) == 2
    @test isa(sol1.alg, ASB07)
    @test sol1.alg.δ == 0.04
    @test sol1.alg.max_order == 5
    @test typeof(sol1.alg.approx_model) == CorrectionHull{EXP.IntervalExpAlg}
    @test sol1.alg.approx_model.order == 10
    @test setrep(sol1) <: Zonotope
    @test setrep(sol1) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # change some of the algorithm-specific parameters
    Calg = CorrectionHull(; order=5,                 # truncation order in Φ
                          exp=EXP.IntervalExpAlg(9)) # truncation order in F
    Ralg = ASB07(; δ=0.01,
                 max_order=7, # maximum zonotope order
                 approx_model=Calg)

    sol2 = solve(ivp; tspan=(0.0, 5.0), alg=Ralg)
    @test dim(sol2) == 2
    @test isa(sol2.alg, ASB07)
    @test sol2.alg.δ == 0.01
    @test sol2.alg.max_order == 7
    @test typeof(sol2.alg.approx_model) == CorrectionHull{EXP.IntervalExpAlg}
    @test sol2.alg.approx_model.order == 5
    @test sol2.alg.approx_model.exp.order == 9
    @test setrep(sol2) <: Zonotope
    @test setrep(sol2) == Zonotope{Float64,Array{Float64,1},Array{Float64,2}}

    # compare recursive option: recursive is more precise
    sol1_rec = solve(ivp; tspan=(0.0, 1.0), alg=ASB07(; δ=0.04, recursive=true))
    sol1_nonrec = solve(ivp; tspan=(0.0, 1.0), alg=ASB07(; δ=0.04, recursive=false))
    @test diameter(set(sol1_rec[end])) < diameter(set(sol1_nonrec[end]))

    # TODO test without IntervalMatrix wrapper
    #A = [-1.0 ± 0.05 -4.0 ± 0.05;
    #     4.0 ± 0.05 -1.0 ± 0.05]
    #ivp = @ivp(x' = A * x, x(0) ∈ X0)
    #sol1 = solve(ivp, tspan=(0.0, 1.0), alg=ASB07(δ=0.04));

    # invariant
    X = BallInf([1.0, 1.0], 100.0)
    ivp = @ivp(x' = A * x, x(0) ∈ X0, x ∈ X)
    sol = solve(ivp; tspan=(0.0, 1.0), alg=alg)
    @test length(sol) == 25
end

@testset "ASB07 algorithm: inhomogeneous case" begin
    # affine ODE: x' = Ax + Bu
    A = [-1.0±0.05 -4.0±0.05;
         4.0±0.05 -1.0±0.05]
    B = hcat([1.0 ± 0.01; 1.0 ± 0.0])
    U = Interval(-0.05, 0.05)
    X0 = BallInf([1.0, 1.0], 0.1)
    ivp = @ivp(x' = A * x + Bu, x(0) ∈ X0, u ∈ U, x ∈ Universe(2))
    alg = ASB07(; δ=0.04)
    sol = solve(ivp; tspan=(0.0, 1.0), alg=alg)
    @test length(sol) == 25

    # invariant
    X = BallInf([1.0, 1.0], 100.0)
    ivp = @ivp(x' = A * x + Bu, x(0) ∈ X0, u ∈ U, x ∈ X)
    sol = solve(ivp; tspan=(0.0, 1.0), alg=alg)
    @test length(sol) == 25
end

@testset "ASB07 with time-triggered hybrid system" begin
    ivp = embrake_pv_1(; ζ=[-1e-8, 1e-7], Tsample=1e-4)
    sol = solve(ivp; alg=ASB07(; δ=1e-7, max_order=3), max_jumps=2)
    @test dim(sol) == 4
    @test flowpipe(sol) isa HybridFlowpipe

    @test tspan(shift(sol[1], 1.0)) == tspan(sol[1]) + 1.0
end
