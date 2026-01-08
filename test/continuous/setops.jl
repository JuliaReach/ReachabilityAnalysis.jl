using Test, ReachabilityAnalysis

@testset "Continuous.Setops" begin
    M = hcat([2.0, 1])
    Z1 = Zonotope([1.0], hcat([1.0]))
    Z2 = Zonotope([3.0, -3], [2.0 0 3; 0 4 5])
    Z3 = ReachabilityAnalysis.linear_map_minkowski_sum(M, Z1, Z2)
    @test Z3 == minkowski_sum(linear_map(M, Z1), Z2)
end
