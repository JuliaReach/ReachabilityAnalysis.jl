@testset "CARLIN algorithm: logistic model" begin
    r = -0.5
    K = 0.8
    F1 = hcat(r)
    F2 = hcat(-r / K)
    X0 = 0.47 .. 0.53
    prob = @ivp(CanonicalQuadraticForm(F1, F2), x(0) ∈ X0)

    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=true))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.2

    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=false))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.01

    # Use the compressed Kronecker form
    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=false, compress=true))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.01
end
