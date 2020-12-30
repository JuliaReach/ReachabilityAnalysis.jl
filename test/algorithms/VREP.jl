@testset "VREP constructor" begin
    N = Float64

    # dimension unknown => VPolytope
    @test setrep(VREP(δ=0.1)) == VPolytope{N, Vector{N}}

    # dimension is two => VPolygon with static arrays
    @test setrep(VREP(δ=0.1, dim=2)) == VPolygon{N, SVector{2, N}}

    # dimension is two => VPolygon with static arrays
    @test setrep(VREP(δ=0.1, dim=2, static=true)) == VPolygon{N, SVector{2, N}}

    # dimension is two, no static
    @test setrep(VREP(δ=0.1, dim=2, static=false)) == VPolygon{N, Vector{N}}

    # dimension is two, use static arrays by default
    @test setrep(VREP(δ=0.1, dim=5)) == VPolytope{N, SVector{5, N}}
end

@testset "VREP static arrays conversions" begin
    N = Float64

    V2 = convert(VPolygon, BallInf(zeros(2), 1.0))
    V3 = convert(VPolytope, BallInf(zeros(3), 1.0))

    # no-op
    @test RA._reconvert(V2, Val(false), Missing) == V2
    @test RA._reconvert(V3, Val(false), Missing) == V3

    # convert to static arrays
    @test RA._reconvert(V2, Val(true), Val(2)) isa VPolygon{Float64, SVector{2, N}}
    @test RA._reconvert(V3, Val(true), Val(3)) isa VPolytope{Float64, SVector{3, N}}
end

@testset "VREP for 2D system" begin
    A  = [0   1.
          -1  0]
    X0 = VPolygon([[1.0, 0.0], [1.2, 0.0], [1.1, 0.2]])
    prob = @ivp(x' = A*x, x(0) ∈ X0)
    sol = solve(prob, tspan=(0.0, 2.0), alg=VREP(δ=1e-2))
    @test setrep(sol) == VPolygon{Float64, Vector{Float64}}

    # statically sized array
    sol = solve(prob, tspan=(0.0, 2.0), alg=VREP(δ=1e-2, static=true, dim=2))
    @test setrep(sol) == VPolygon{Float64, SVector{2, Float64}}
end

@testset "VREP for 5D system" begin
    A  = [1. 0  0 -5
          4 -2  4 -3
         -4  0  0  1
          5 -2  2  3.]

    X0_12 = VPolygon([[1.0, 0.0], [1.2, 0.0], [1.1, 0.2]])
    X0_34 = BallInf(zeros(2), 0.2)
    X0 = VPolytope(vertices_list(X0_12 × X0_34))

    # default Polyhedra backend works, but shows warning messages of this sort:
    # glp_simplex: unable to recover undefined or non-optimal solution
    #prob = @ivp(x' = A*x, x(0) ∈ X0)
    #sol = solve(prob, tspan=(0.0, 1.0), alg=VREP(δ=1e-3, static=true, dim=4))

    # change the backend
    prob = @ivp(x' = A*x, x(0) ∈ X0)
    sol = solve(prob, tspan=(0.0, 1.0), alg=VREP(δ=1e-3, static=true, dim=4, backend=CDDLib.Library()))
end
