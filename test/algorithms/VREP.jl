using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
import CDDLib
using StaticArrays: SVector

@testset "VREP struct" begin
    N = Float64

    # full constructor
    alg = VREP(0.01, nothing, nothing, nothing)

    # constructor with default values
    alg = VREP(; δ=0.01)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == N
    @test setrep(alg) == VPolytope{N,Vector{N}}
    @test rsetrep(alg) == ReachSet{N,VPolytope{N,Vector{N}}}

    # dimension is two => VPolygon with static arrays by default
    @test setrep(VREP(; δ=0.1, dim=2)) == VPolygon{N,SVector{2,N}}
    @test setrep(VREP(; δ=0.1, dim=2, static=false)) == VPolygon{N,Vector{N}}

    # dimension is not two => VPolytope with static arrays by default
    @test setrep(VREP(; δ=0.1, dim=5)) == VPolytope{N,SVector{5,N}}
    @test_broken setrep(VREP(; δ=0.1, dim=5, static=false)) == VPolytope{N,Vector{N}}  # TODO fix this

    # `static` option requires `dim`
    @test_throws ErrorException setrep(VREP(; δ=0.1, static=true))
end

@testset "VREP static arrays conversions" begin
    N = Float64

    V2 = convert(VPolygon, BallInf(zeros(2), 1.0))
    V3 = convert(VPolytope, BallInf(zeros(3), 1.0))

    # no-op
    @test RA._reconvert(V2, Val(false), Missing) == V2
    @test RA._reconvert(V3, Val(false), Missing) == V3

    # convert to static arrays
    @test RA._reconvert(V2, Val(true), Val(2)) isa VPolygon{Float64,SVector{2,N}}
    @test RA._reconvert(V3, Val(true), Val(3)) isa VPolytope{Float64,SVector{3,N}}
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
    ivp_norm = RA._normalize(ivp)
    ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
    NSTEPS = 500
    fp_d = RA.post(alg, ivp_discr, NSTEPS)

    # statically sized array
    sol = solve(ivp; tspan=(0.0, 2.0), alg=VREP(; δ=1e-2, static=true, dim=2))
    @test setrep(sol) == VPolygon{Float64,SVector{2,Float64}}

    # inhomogeneous system
    B = [1.0 0; 0 1]
    U = BallInf(zeros(2), 0.1)
    X = Universe(2)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)
    @test_broken discretize(normalize(ivpI), alg.δ, alg.approx_model)  # TODO implement this
    @test_broken solve(ivpI; tspan=(0.0, 2.0), alg=alg) isa Flowpipe  # TODO implement this
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
