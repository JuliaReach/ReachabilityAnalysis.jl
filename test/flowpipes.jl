@testset "Abstract flowpipe interface" begin
    prob, dt = harmonic_oscillator()
    sol = solve(prob, tspan=dt)
    fp = flowpipe(sol)

    FT = typeof(fp)
    @test basetype(FT) == basetype(fp) == Flowpipe

    d = ones(2)
    @test ρ(d, fp) < 2.0

end


@testset "Flowpipe constructors" begin

    # uninitialized array constructor
    # if the set type is passed, then a ReachSet is used by default
    fp = Flowpipe(undef, Singleton{Float64, Vector{Float64}}, 3)
    @test length(fp) == 3
    @test_throws UndefRefError dim(fp) # fp is not empty but its dimension is not defined
    @test setrep(fp) == Singleton{Float64, Vector{Float64}}

    # the reach-set type can be passed directly
    fpr = Flowpipe(undef, ReachSet{Float64, Singleton{Float64, Vector{Float64}}}, 3)
    @test length(fpr) == 3
    @test typeof(fpr) == typeof(fp)
    @test setrep(fpr) == Singleton{Float64, Vector{Float64}}

    # try another reach-set type
    RSP = SparseReachSet{Float64, Singleton{Float64, Vector{Float64}}, 2}
    fp = Flowpipe(undef, RSP, 3)
    @test length(fp) == 3
    @test dim(fp) == 2 # in this set type, the dimension is known
    @test setrep(fp) == Singleton{Float64, Vector{Float64}}
end

@testset "Hybrid flowpipe time evaluation" begin
    X = Interval(0 .. 1)
    δ = 0.1
    F1 = Flowpipe([ReachSet(X, (0 .. δ) + k*δ) for k in 0:10])
    F2 = Flowpipe([ReachSet(X, (0 .. δ) + k*δ) for k in 9:20])
    F2′ = Flowpipe([ReachSet(X, (0 .. δ) + k*δ*0.9) for k in 9:20])
    H = HybridFlowpipe([F1, F2])
    H′ = HybridFlowpipe([F1, F2′])

    times = 0:0.01:1.9
    values = [[x] for x in (@. 1 - exp(-times))]
    @test all(vi ∈ H(ti) for (ti, vi) in zip(times, values))
    @test all(vi ∈ H′(ti) for (ti, vi) in zip(times, values))
end

@testset "Flowpipe clustering" begin
    X = Interval(0 .. 1)
    δ = 0.1
    F1 = Flowpipe([ReachSet(X, (0 .. δ) + k*δ) for k in 0:10])

    # FIXME requires LazySets#2157
    #N = eltype(X)
    #U1 = cluster(F1, 8:10, UnionClustering())
    #@test U1 isa Vector{ReachSet{N, UnionSetArray{N,Interval{N, IA.Interval{N}}}}}
    #@test set(first(U1)) == UnionSetArray([set(Fi) for Fi in F1[8:10]])
end
