@testset "Abstract flowpipe interface" begin
    prob, dt = harmonic_oscillator()
    sol = solve(prob; tspan=dt)
    fp = flowpipe(sol)

    # base, numeric, set and reach-set types
    FT = typeof(fp)
    @test basetype(FT) == basetype(fp) == Flowpipe
    N = ReachabilityAnalysis.numtype(fp)
    @test N == Float64
    ZT = Zonotope{N,Vector{N},Matrix{N}}
    @test setrep(FT) == ZT
    @test rsetrep(fp) == ReachSet{N,ZT}

    # ambient dimension
    @test dim(fp) == 2

    # indexing
    @test first(fp) == fp[1] === fp.Xk[1]
    @test last(fp) == fp[end] === fp.Xk[end]
    @test length(fp) == length(fp.Xk) == numrsets(fp)
    @test set(fp, 3) === set(fp[3])
    @test set(fp, 3:5) isa UnionSetArray
    @test set(fp) isa UnionSetArray
    @test array(fp) === fp.Xk

    # time domain interface
    @test tstart(fp) ≈ 0.0 && dt[1] ≈ 0.0
    @test tend(fp) ≈ 20.0 && dt[2] ≈ 20.0
    @test RA._isapprox(tspan(fp), 0 .. 20)
    @test vars(fp) == (1, 2)
end

@testset "Flowpipe set operations interface" begin
    prob, dt = harmonic_oscillator()
    sol = solve(prob; tspan=dt)
    fp = flowpipe(sol)

    # support function and support vector
    d = ones(2)
    @test ρ(d, fp) < 2.0
    @test σ(d, fp) ≈ [0.7502062459270877, 0.6989729206916551]

    # disjointness
    H = HalfSpace(d, 2.0)
    Hc = complement(H)

    # inclusion

    # membership

end

@testset begin
    "Flowpipe inclusion"
    # Projectile problem 
    A = [0.0 0.5 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.7; 0.0 0.0 0.0 0.0]
    X0 = Singleton([0.0, 5.0, 100.0, 0])
    U = Singleton([0.0, 0.0, 0.0, -9.81])
    prob = @ivp(x' = A * x + Matrix(1.0I, 4, 4) * u, x(0) ∈ X0, u ∈ U, x ∈ Universe(4))
    cons = LinearConstraint([24.0, 0.0, 1, 0], 375.0)
    alg = GLGM06(; δ=1e-2, approx_model=Forward())
    sol = solve(prob; T=20.0, alg=alg)

    # equivalent ways
    @test all(set(R) ⊆ cons for R in sol)
    @test all(R ⊆ cons for R in sol)
    @test sol ⊆ cons
end

@testset "Flowpipe constructors" begin
    # uninitialized array constructor
    # if the set type is passed, then a ReachSet is used by default
    fp = Flowpipe(undef, Singleton{Float64,Vector{Float64}}, 3)
    @test length(fp) == 3
    @test_throws UndefRefError dim(fp) # fp is not empty but its dimension is not defined
    @test setrep(fp) == Singleton{Float64,Vector{Float64}}

    # the reach-set type can be passed directly
    fpr = Flowpipe(undef, ReachSet{Float64,Singleton{Float64,Vector{Float64}}}, 3)
    @test length(fpr) == 3
    @test typeof(fpr) == typeof(fp)
    @test setrep(fpr) == Singleton{Float64,Vector{Float64}}

    # try another reach-set type
    RSP = SparseReachSet{Float64,Singleton{Float64,Vector{Float64}},2}
    fp = Flowpipe(undef, RSP, 3)
    @test length(fp) == 3
    @test dim(fp) == 2 # in this set type, the dimension is known
    @test setrep(fp) == Singleton{Float64,Vector{Float64}}
end

@testset "Hybrid flowpipe time evaluation" begin
    X = Interval(0 .. 1)
    δ = 0.1
    F1 = Flowpipe([ReachSet(X, (0 .. δ) + k * δ) for k in 0:10])
    F2 = Flowpipe([ReachSet(X, (0 .. δ) + k * δ) for k in 9:20])
    F2′ = Flowpipe([ReachSet(X, (0 .. δ) + k * δ * 0.9) for k in 9:20])
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
    F1 = Flowpipe([ReachSet(X, (0 .. δ) + k * δ) for k in 0:10])

    #= FIXME requires LazySets#2157
    N = eltype(X)
    U1 = cluster(F1, 8:10, UnionClustering())
    @test U1 isa Vector{ReachSet{N, UnionSetArray{N,Interval{N, IA.Interval{N}}}}}
    @test set(first(U1)) == UnionSetArray([set(Fi) for Fi in F1[8:10]])

    # Test ReachabilityAnalysis#347
    prob, dt = vanderpol()
    sol = solve(prob, tspan=(0.0, 3.0));
    fp = flowpipe(sol)
    cl = cluster(fp, [1, 2, 3], UnionClustering());
    @test length(cl) == 1 && typeof(set(first(cl))) <: UnionSetArray
    =#
end
