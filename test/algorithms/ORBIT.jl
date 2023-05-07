using LazySets: _isapprox

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
    @test _isapprox(tspan(sol), interval(0, 3δ))
    x = sol[end]

    @test _isapprox(tspan(x), interval(3δ))
    @test element(set(x)) ≈ [cos(3δ), -sin(3δ)]

    # test time span sequence
    @test all(_isapprox(tspan(sol[i]), interval(δ * (i - 1))) for i in eachindex(sol))

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
