using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA

@testset "BFFPSV18 struct" begin
    # full constructor (requires to manually define the type parameter for the set representation)
    alg = BFFPSV18{Float64,Int,Int,Int,Int,Int,Int}(0.01, 0, 0, 0, 0, 0, false, false, false, false)

    # constructor with default values
    alg = BFFPSV18(; δ=0.01, dim=1)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test setrep(alg) == Interval{Float64}
    @test rsetrep(alg) == SparseReachSet{Float64,CartesianProductArray{Float64,setrep(alg)},1}

    # higher-dimensional blocks
    alg = BFFPSV18(; δ=0.01, dim=2)
    @test setrep(alg) == Interval{Float64}
    @test rsetrep(alg) == SparseReachSet{Float64,CartesianProductArray{Float64,setrep(alg)},2}
    alg = BFFPSV18(; δ=0.01, partition=[1:2, 3:4])
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    @test rsetrep(alg) == SparseReachSet{Float64,CartesianProductArray{Float64,setrep(alg)},4}

    # invalid construction: requires `dim` or `partition`
    @test_throws ArgumentError BFFPSV18(; δ=0.01)
end

@testset "BFFPSV18 algorithm: type constructors" begin
    A = rand(5, 5)
    X0 = rand(Hyperrectangle; dim=5)
    B = rand(5, 1)
    X = Universe(5)
    U = rand(Hyperrectangle; dim=1)
    ivp = @ivp(x' = A * x, x(0) ∈ X0)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, x ∈ X, u ∈ U)
    δ = 1e-3
    T = 2 * δ

    # interval (1D) decomposition, one variable
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=[2], dim=5)
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test setrep(alg) == Interval{Float64}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # interval (1D) decomposition, one variable, without explicitly passing Interval
    alg = BFFPSV18(; δ=δ, vars=[2], dim=5)
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test setrep(alg) == Interval{Float64}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # interval (1D) decomposition, two variables
    alg = BFFPSV18(; δ=δ, vars=[1, 4], dim=5)
    @test alg.vars == [1, 4]
    @test alg.block_indices == [1, 4]
    @test alg.row_blocks == [[1], [4]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test setrep(alg) == Interval{Float64}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangle decomposition
    alg = BFFPSV18(; δ=δ, setrep=Hyperrectangle, vars=[3], partition=[[1], [2, 3, 4, 5]])
    @test alg.vars == [3]
    @test alg.block_indices == [2] # variable 3 is in the second block
    @test alg.row_blocks == [[2, 3, 4, 5]]
    @test alg.column_blocks == [[1], [2, 3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangle decomposition, without explicitly passing the set type
    alg = BFFPSV18(; δ=δ, vars=[3], partition=[[1], [2, 3, 4, 5]])
    @test alg.vars == [3]
    @test alg.block_indices == [2] # variable 3 is in the second block
    @test alg.row_blocks == [[2, 3, 4, 5]]
    @test alg.column_blocks == [[1], [2, 3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangular decomposition with different block sizes and a single target variable
    alg = BFFPSV18(; δ=δ, setrep=Hyperrectangle, vars=[3], partition=[1:2, 3:3, 4:5])
    @test alg.vars == [3]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [3:3]
    @test alg.column_blocks == [1:2, 3:3, 4:5]
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangular decomposition with different block sizes, with variables belonging to
    # different blocks
    alg = BFFPSV18(; δ=δ, setrep=Hyperrectangle, vars=[3, 4, 5], partition=[1:2, 3:3, 4:5])
    @test alg.vars == [3, 4, 5]
    @test alg.block_indices == [2, 3]
    @test alg.row_blocks == [3:3, 4:5]
    @test alg.column_blocks == [1:2, 3:3, 4:5]
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangle decomposition with different variables in the same partition
    alg = BFFPSV18(; δ=δ, vars=[1, 2], partition=[1:2, 3:5])
    @test alg.vars == [1, 2]
    @test alg.block_indices == [1]
    @test alg.row_blocks == [[1, 2]]
    @test alg.column_blocks == [[1, 2], [3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # hyperrectangle decomposition with different variables in different partitions
    alg = BFFPSV18(; δ=δ, vars=[1, 3], partition=[1:2, 3:5])
    @test alg.vars == [1, 3]
    @test alg.block_indices == [1, 2]
    @test alg.row_blocks == [[1, 2], [3, 4, 5]]
    @test alg.column_blocks == [[1, 2], [3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64,Vector{Float64},Vector{Float64}}
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # sparse option, interval (1D) decomposition
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=[2], dim=5, sparse=true)
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test alg.sparse == true
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # sparse option, generic (nD) decomposition
    alg = BFFPSV18(; δ=δ, setrep=Hyperrectangle, vars=[2], dim=5, sparse=true)
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test alg.sparse == true
    sol = solve(ivp; T=T, alg=alg)
    sol = solve(ivpI; T=T, alg=alg)

    # sparse option, no view
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=[2], dim=5, sparse=true, view=false)
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4], [5]]
    @test alg.sparse == true
    @test alg.view == false
    @test_broken solve(ivp; T=T, alg=alg)  # TODO missing method
    sol = solve(ivpI; T=T, alg=alg)

    # invariant
    X = BallInf(zeros(5), 100.0)
    ivp = @ivp(x' = A * x, x(0) ∈ X0, x ∈ X)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, x ∈ X, u ∈ U)
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=1:5, dim=5)
    @test_broken solve(ivp; T=T, alg=alg)  # TODO missing method
    sol = solve(ivpI; T=T, alg=alg)

    # invariant incompatible with selected variables
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=[2], dim=5)
    @test_broken solve(ivp; T=T, alg=alg)  # TODO missing method
    @test_throws DimensionMismatch solve(ivpI; T=T, alg=alg)

    # invariant compatible with selected variable
    X = BallInf(zeros(1), 100.0)
    ivp = @ivp(x' = A * x, x(0) ∈ X0, x ∈ X)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, x ∈ X, u ∈ U)
    alg = BFFPSV18(; δ=δ, setrep=Interval, vars=[2], dim=5)
    @test_broken solve(ivp; T=T, alg=alg)  # TODO missing method
    sol = solve(ivpI; T=T, alg=alg)
end
