import Polyhedra
import CDDLib
using StaticArrays: SVector
using ReachabilityAnalysis: _reconvert

@testset "VREP constructor" begin
    N = Float64

    # dimension unknown => VPolytope
    @test setrep(VREP(; δ=0.1)) == VPolytope{N,Vector{N}}

    # dimension is two => VPolygon with static arrays
    @test setrep(VREP(; δ=0.1, dim=2)) == VPolygon{N,SVector{2,N}}

    # dimension is two => VPolygon with static arrays
    @test setrep(VREP(; δ=0.1, dim=2, static=true)) == VPolygon{N,SVector{2,N}}

    # dimension is two, no static
    @test setrep(VREP(; δ=0.1, dim=2, static=false)) == VPolygon{N,Vector{N}}

    # dimension is two, use static arrays by default
    @test setrep(VREP(; δ=0.1, dim=5)) == VPolytope{N,SVector{5,N}}
end

@testset "VREP static arrays conversions" begin
    N = Float64

    V2 = convert(VPolygon, BallInf(zeros(2), 1.0))
    V3 = convert(VPolytope, BallInf(zeros(3), 1.0))

    # no-op
    @test _reconvert(V2, Val(false), Missing) == V2
    @test _reconvert(V3, Val(false), Missing) == V3

    # convert to static arrays
    @test _reconvert(V2, Val(true), Val(2)) isa VPolygon{Float64,SVector{2,N}}
    @test _reconvert(V3, Val(true), Val(3)) isa VPolytope{Float64,SVector{3,N}}
end

@testset "VREP for 2D system" begin
    A = [0 1.0
         -1 0]
    X0 = VPolygon([[1.0, 0.0], [1.2, 0.0], [1.1, 0.2]])
    ivp = @ivp(x' = A * x, x(0) ∈ X0)
    alg = VREP(; δ=1e-2)
    # continuous algorithm
    sol = solve(ivp; tspan=(0.0, 2.0), alg=alg)
    @test setrep(sol) == VPolygon{Float64,Vector{Float64}}
    # discrete algorithm
    ivp_norm = ReachabilityAnalysis._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

    # statically sized array
    sol = solve(ivp; tspan=(0.0, 2.0), alg=VREP(; δ=1e-2, static=true, dim=2))
    @test setrep(sol) == VPolygon{Float64,SVector{2,Float64}}
end

@testset "VREP for 5D system" begin
    A = [1.0 0 0 -5
         4 -2 4 -3
         -4 0 0 1
         5 -2 2 3.0]

    X0_12 = VPolygon([[1.0, 0.0], [1.2, 0.0], [1.1, 0.2]])
    X0_34 = BallInf(zeros(2), 0.2)
    X0 = VPolytope(vertices_list(X0_12 × X0_34))

    # default Polyhedra backend works, but shows warning messages of this sort:
    # glp_simplex: unable to recover undefined or non-optimal solution
    #ivp = @ivp(x' = A*x, x(0) ∈ X0)
    #sol = solve(ivp, tspan=(0.0, 1.0), alg=VREP(δ=1e-3, static=true, dim=4))

    # change the backend
    ivp = @ivp(x' = A * x, x(0) ∈ X0)
    sol = solve(ivp; tspan=(0.0, 1.0),
                alg=VREP(; δ=1e-3, static=true, dim=4, backend=CDDLib.Library()))
end

@testset "VREP for sdof with distinct approximations" begin
    ivp, _ = harmonic_oscillator()
    tmax = 2.0

    a = solve(ivp; tspan=(0.0, tmax), alg=VREP(; δ=0.1, approx_model=Forward(; setops=:concrete)))
    b = solve(ivp; tspan=(0.0, tmax),
              alg=VREP(; δ=0.1, approx_model=StepIntersect(; setops=:concrete)))
    c = solve(ivp; tspan=(0.0, tmax), alg=VREP(; δ=0.1, approx_model=CorrectionHull(; exp=:base)))
end
