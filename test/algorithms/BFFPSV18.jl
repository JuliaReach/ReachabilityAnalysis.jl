@testset "BFFPSV18 algorithm: type constructors" begin

    # interval (1D) decomposition, one variable
    alg = BFFPSV18(δ=1e-3, setrep=Interval, vars=[2], dim=4);
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4]]
    @test setrep(alg) == Interval{Float64, IA.Interval{Float64}}

    # interval (1D) decomposition, one variable, without explicitly passing Interval
    alg = BFFPSV18(δ=1e-3, vars=[2], dim=4);
    @test alg.vars == [2]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[2]]
    @test alg.column_blocks == [[1], [2], [3], [4]]
    @test setrep(alg) == Interval{Float64, IA.Interval{Float64}}

    # interval (1D) decomposition, two variables
    alg = BFFPSV18(δ=1e-3, vars=[1, 4], dim=4);
    @test alg.vars == [1, 4]
    @test alg.block_indices == [1, 4]
    @test alg.row_blocks == [[1], [4]]
    @test alg.column_blocks == [[1], [2], [3], [4]]
    @test setrep(alg) == Interval{Float64, IA.Interval{Float64}}

    # hyperrectangle decomposition
    alg = BFFPSV18(δ=1e-3, setrep=Hyperrectangle, vars=[3], partition=[[1], [2, 3, 4]]);
    @test alg.vars == [3]
    @test alg.block_indices == [2] # variable 3 is in the second block
    @test alg.row_blocks == [[2, 3, 4]]
    @test alg.column_blocks == [[1], [2, 3, 4]]
    @test setrep(alg) == Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}

    # hyperrectangle decomposition, without explicitly passing the set type
    alg = BFFPSV18(δ=1e-3, vars=[3], partition=[[1], [2, 3, 4]]);
    @test alg.vars == [3]
    @test alg.block_indices == [2] # variable 3 is in the second block
    @test alg.row_blocks == [[2, 3, 4]]
    @test alg.column_blocks == [[1], [2, 3, 4]]
    @test setrep(alg) == Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}

    # hyperrectangular decomposition with different block sizes and a single target variable
    alg = BFFPSV18(δ=1e-3, setrep=Hyperrectangle, vars=[25], partition = [1:24, [25], 26:48])
    @test alg.vars == [25]
    @test alg.block_indices == [2]
    @test alg.row_blocks == [[25]]
    @test alg.column_blocks == [1:24, [25], 26:48]

    # hyperrectangular decomposition with different block sizes, with varables belonging to different blocks
    alg = BFFPSV18(δ=1e-3, setrep=Hyperrectangle, vars=[25, 37, 40], partition = [1:24, [25], 26:48])
    @test alg.vars == [25, 37, 40]
    @test alg.block_indices == [2, 3]
    @test alg.row_blocks == [[25], 26:48]
    @test alg.column_blocks == [1:24, [25], 26:48]

    # hyperrectangle decomposition with different variables in the same partition
    alg = BFFPSV18(δ=1e-3, vars=[1, 2], partition=[1:2, 3:5]);
    @test alg.vars == [1, 2]
    @test alg.block_indices == [1]
    @test alg.row_blocks == [[1, 2]]
    @test alg.column_blocks == [[1, 2], [3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}

    # hyperrectangle decomposition with different variables in different partitions
    alg = BFFPSV18(δ=1e-3, vars=[1, 3], partition=[1:2, 3:5]);
    @test alg.vars == [1, 3]
    @test alg.block_indices == [1, 2]
    @test alg.row_blocks == [[1, 2], [3, 4, 5]]
    @test alg.column_blocks == [[1, 2], [3, 4, 5]]
    @test setrep(alg) == Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}
end
