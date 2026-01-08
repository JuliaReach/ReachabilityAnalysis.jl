@testset "discretization: CorrectionHull" begin
    A = IntervalMatrix([interval(-1.05, -0.95) interval(-4.05, -3.95);
                        interval(3.95, 4.05) interval(-1.05, -0.95)])
    B = hcat([interval(0.99, 1.01); interval(1)])
    X0 = BallInf([1.0, 1.0], 0.1)
    δ = 0.01
    alg = CorrectionHull()

    # origin is contained in U
    U = Interval(-0.05, 0.05)
    ivp = @ivp(x' = A * x + Bu, x(0) ∈ X0, u ∈ U, x ∈ Universe(2))
    ivp = normalize(ivp)
    discretize(ivp, δ, alg)

    # origin is not contained in U
    U = Interval(0.05, 0.1)
    ivp = @ivp(x' = A * x + Bu, x(0) ∈ X0, u ∈ U, x ∈ Universe(2))
    ivp = normalize(ivp)
    discretize(ivp, δ, alg)
end
