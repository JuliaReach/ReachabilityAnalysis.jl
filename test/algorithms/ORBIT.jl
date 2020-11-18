@testset "ORBIT algorithm" begin
    δ = 1e-2

    # harmonic oscillator x'' = -x, s.t. x(0) = 1, x'(0) = 0
    prob, dt = harmonic_oscillator()
    sol = solve(prob, tspan=(0.0, 20.0), alg=ORBIT(δ=δ));

    # check result
    @test sum((set(sol[k]).element[1] - cos(δ * (k-1)))^2 for k in 1:length(sol)) < eps(Float64)
    @test sum((set(sol[k]).element[2] + sin(δ * (k-1)))^2 for k in 1:length(sol)) < eps(Float64)

    # forced oscillator Mx'' + Kx = R, s.t. x(0) = [0, 0], x'(0) = [0, 0]
    prob, dt = forced_oscillator()
    sol = solve(prob, tspan=(0.0, 20.0), alg=ORBIT(δ=δ))

    # check result
    U = forced_oscillator_solution()
    U1 = k -> U((k-1) * δ)[1]
    U2 = k -> U((k-1) * δ)[2]
    @test sum((set(sol[k]).element[1] - U1(k))^2 for k in 1:length(sol)) < eps(Float64)
    @test sum((set(sol[k]).element[2] - U2(k))^2 for k in 1:length(sol)) < eps(Float64)
end