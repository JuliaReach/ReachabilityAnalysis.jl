@testset "Flowpipe interface" begin
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    sol = solve(p, tspan=(0.0, 1.0), GLGM06(δ=0.1))
    F = flowpipe(sol)

    # iteration
    @test set(F, 1) == set(F[1])
    @test set(F, length(F)) == set(F[end])

    # time span methods
    @test _isapprox(tspan(sol), 0.0 .. 1.0)
    @test tstart(F) == 0.0
    @test tend(F) ≈ 1.0
    @test tstart(F[1]) == 0.0
    @test tend(F[1]) ≈ 0.1
    @test tstart(F[end]) ≈ 0.9
    @test tend(F[end]) ≈ 1.0

    # set representation
    #@test setrep(F) isa Zonotope
    @test setrep(F) == Zonotope{Float64, Vector{Float64}, Matrix{Float64}}

    # callable behavior
    @test F(0..1) == F[1:end]
    @test F(0.05) == F[1]
    @test F(1.0) == F[end]
    @test F(0.05 .. 0.15) == F[1:2]
    @test_throws ArgumentError F(0.05 .. 1.05) # time interval not contained in flowpipe
end

@testset "Solution interface: initial states" begin
    # interval initial condition
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    solve(p, T=1.0)

    # interval initial condition
    p = @ivp(x' = -x, x(0) ∈ Interval(0, 1))
    solve(p, T=1.0)

    # deterministic initial condition, scalar for one-dimensional matrix (TODO)
    #p = InitialValueProblem(@system(x' = -x), 0.5)
    #solve(p, T=1.0)

    # deterministic initial condition, vector
    #p = InitialValueProblem(@system(x' = -x), [0.5])
    #solve(p, T=1.0)
end

@testset "Solution interface: time span" begin
    p = @ivp(x' = -x, x(0) ∈ 0..1)
    Δt = 0.0 .. 2.0

    # only time horizon given
    sol = solve(p, T=2.0)
    @test _isapprox(tspan(sol), Δt)

    # both time horizon and time span given (TODO)
    #@test_throws ArgumentError solve(p, T=1.0, timespan=(0.0, 2.0))

    # time span given as a tuple
    sol = solve(p, tspan=(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as an interval
    sol = solve(p, tspan=0.0 .. 2.0)
    @test _isapprox(tspan(sol), Δt)

    # time span given as a LazySets interval
    sol = solve(p, tspan=Interval(0.0, 2.0))
    @test _isapprox(tspan(sol), Δt)

    # time span given as a vector
    sol = solve(p, tspan=[0.0, 2.0])
    @test _isapprox(tspan(sol), Δt)
end
