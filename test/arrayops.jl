using ReachabilityAnalysis: logarithmic_norm

@testset "Logarithmic norm" begin

    # Example 1
    A = [1 2; 3 -4]

    @test logarithmic_norm(A, 1) == 4
    @test logarithmic_norm(A, Inf) == 3

    # Example 2 (see Example 2 in "Explicit Error Bounds for Carleman Linearization.
    # Marcelo Forets and Amaury Pouly (2018)").
    A = [0 1; -1 -2]
    Aop2 = kron(A, A)

    @test logarithmic_norm(A, Inf) == 1
    @test logarithmic_norm(Aop2, Inf) == 9;

    @test logarithmic_norm(A, 1) == 1
    @test logarithmic_norm(Aop2, 1) == 9

    @test logarithmic_norm(A, 2) == 0
    @test logarithmic_norm(Aop2, 2) ≈ 4.23606797749979
end
