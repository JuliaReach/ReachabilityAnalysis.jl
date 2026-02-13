using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
using ReachabilityAnalysis: _isapprox

@testset "ORBIT struct" begin
    # full constructor (requires to manually define the type parameter for the set representation)
    alg = ORBIT{Float64,Int,Int}(0.01, 0)

    # constructor with default values
    alg = ORBIT(; δ=0.01)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test setrep(alg) == Singleton{Float64,Vector{Float64}}
    @test rsetrep(alg) == ReachSet{Float64,setrep(alg)}
end

@testset "ORBIT algorithm" begin
    # --------------------------------------------------------------------------
    # harmonic oscillator x'' = -x, s.t. x(0) = 1, x'(0) = 0
    # --------------------------------------------------------------------------
    prob, dt = harmonic_oscillator()
    δ = 1e-2
    sol = solve(prob; tspan=(0.0, 20.0), alg=ORBIT(; δ=δ))

    # check result
    @test sum((set(sol[k]).element[1] - cos(δ * (k - 1)))^2 for k in eachindex(sol)) < eps(Float64)
    @test sum((set(sol[k]).element[2] + sin(δ * (k - 1)))^2 for k in eachindex(sol)) < eps(Float64)

    # test number of steps
    sol = solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ))
    @test _isapprox(tspan(sol), TimeInterval(IA.interval(0, 3δ)))
    x = sol[end]

    @test _isapprox(tspan(x), TimeInterval(IA.interval(3δ)))
    @test element(set(x)) ≈ [cos(3δ), -sin(3δ)]

    # test time span sequence
    @test all(_isapprox(tspan(sol[i]), TimeInterval(IA.interval(δ * (i - 1))))
              for i in eachindex(sol))

    # `homogenize` option
    @test_broken solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ), homogenize=true) isa
                 RA.ReachSolution  # TODO fix this

    # start from a ZeroSet
    prob0 = IVP(prob.s, ZeroSet(2))
    @test_broken solve(prob0; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ)) isa RA.ReachSolution  # TODO fix this

    # Krylov
    if isdefined(@__MODULE__, :ExponentialUtilities)
        A = prob.s.A
        res = RA.ORBITModule._orbit_krylov!(A, [1.0, 1.0], 2)
        @test length(res) == 2
    end

    # --------------------------------------------------------------------------
    # inhomogeneous problem
    # --------------------------------------------------------------------------
    A = prob.s.A
    X0 = prob.x0
    B = [2.0 0; 0 1]
    X = Universe(2)

    # AbstractInput
    U = ConstantInput(Singleton([0.1, 0.1]))
    prob = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)
    sol = solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ)) isa RA.ReachSolution

    # start from a ZeroSet and U is a Singleton
    prob = @ivp(x' = A * x + u, x(0) ∈ ZeroSet(2), u ∈ U, x ∈ X)
    sol = solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ))

    # nondeterministic inputs
    U = BallInf([0.1, 0.1], 0.0)
    prob = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)
    @test_broken solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ)) isa RA.ReachSolution  # TODO fix this

    # invariant
    X = BallInf(zeros(2), 100.0)
    prob = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)
    @test_broken solve(prob; tspan=(0.0, 3δ), alg=ORBIT(; δ=δ)) isa RA.ReachSolution  # TODO implement this

    # --------------------------------------------------------------------------
    # forced oscillator Mx'' + Kx = R, s.t. x(0) = [0, 0], x'(0) = [0, 0]
    # --------------------------------------------------------------------------
    prob, dt = forced_oscillator()
    δ = 1e-2
    sol = solve(prob; tspan=(0.0, 20.0), alg=ORBIT(; δ=δ))

    # check result
    U = forced_oscillator_solution()
    U1 = k -> U((k - 1) * δ)[1]
    U2 = k -> U((k - 1) * δ)[2]
    @test sum((set(sol[k]).element[1] - U1(k))^2 for k in eachindex(sol)) < eps(Float64)
    @test sum((set(sol[k]).element[2] - U2(k))^2 for k in eachindex(sol)) < eps(Float64)
end
