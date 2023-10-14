# using SparseArrays
# using Expokit
# using ReachabilityAnalysis: _exp, BaseExpAlg, LazyExpAlg, IntervalExpAlg, PadeExpAlg
using ReachabilityAnalysis.ExponentiationModule

@testset "Exponentiation" begin
    A = [10.0 0; 0 20]
    δ = 0.1
    expAδ = [exp(1) 0; 0 exp(2)]

    # default algorithm
    @test _exp(A, δ) == expAδ

    # base algorithm
    @test _exp(A, δ, BaseExpAlg()) == expAδ

    # Krylov algorithm
    @test _exp(sparse(A), δ, LazyExpAlg()) == expAδ

    # interval matrix
    @test all((_exp(A, δ, IntervalExpAlg()) .- IntervalMatrix(expAδ)) .<= 1e-4)

    # Padé approximation
    @test_throws ArgumentError _exp(A, δ, PadeExpAlg())
    @test _exp(sparse(A), δ, PadeExpAlg()) ≈ expAδ
end
