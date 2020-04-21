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
