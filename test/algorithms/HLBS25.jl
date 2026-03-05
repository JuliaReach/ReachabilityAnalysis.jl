@testset "HLBS25 algorithm: homogeneous case" begin
    for N in (Float32, Float64)

        # uncertain system matrix
        AS = MatrixZonotope(N[-1 -5; 1 -1], [N[0 0; 0.1 0]], [6])

        # initial set: box -> zonotope -> sparse polynomial zonotope
        X0 = Zonotope(N[1.0, 1.0], N[0.2 0; 0 0.2])
        X0 = convert(SparsePolynomialZonotope, X0)

        # problem definition
        prob = @ivp(x' = A * x, x(0) ∈ X0, A ∈ AS)

        T = 2 * N(π)
        δ = T / 400

        # recursive false
        alg = HLBS25(; δ=δ,
                     approx_model=CorrectionHullMatrixZonotope(),
                     taylor_order=5,
                     max_order=5,
                     reduction_method=LazySets.GIR05(),
                     recursive=false)
        sol1 = solve(prob, alg; T=T)

        @test_broken X0 ⊆ set(sol[1]) #TODO add inclusion check for SPZ

        @test isa(sol1.alg, HLBS25)
        @test sol1.alg.δ == δ
        @test sol1.alg.max_order == 5
        @test sol1.alg.max_order_zono == 5
        @test sol1.alg.taylor_order == 5
        @test sol1.alg.approx_model isa CorrectionHullMatrixZonotope
        @test setrep(sol1) == setrep(alg) ==
              SparsePolynomialZonotope{N,Vector{N},Matrix{N},Matrix{N},Matrix{Int},Vector{Int}}
        @test rsetrep(alg) == ReachSet{N,setrep(alg)}

        # recursive true
        alg2 = HLBS25(; δ=δ,
                      approx_model=CorrectionHullMatrixZonotope(; recursive=true),
                      taylor_order=4,
                      max_order=4,
                      reduction_method=LazySets.GIR05(),
                      recursive=true)
        sol2 = solve(prob, alg2; T=T)

        @test isa(sol2.alg, HLBS25)
        @test sol2.alg.δ == δ
        @test sol2.alg.max_order == 4
        @test sol2.alg.max_order_zono == 4
        @test sol2.alg.taylor_order == 4
        @test setrep(sol2) == setrep(alg) ==
              SparsePolynomialZonotope{N,Vector{N},Matrix{N},Matrix{N},Matrix{Int},Vector{Int}}
        @test rsetrep(alg) == ReachSet{N,setrep(alg)}
    end
end

@testset "HLBS25 algorithm: inhomogeneous case" begin
    for N in (Float32, Float64)
        # uncertain system and input matrices
        AS = MatrixZonotope(N[-1 -5; 1 -1], [N[0 0; 0.1 0]], [6])
        BS = MatrixZonotope(N[1 0; 0 1], [N[0.05 0; 0 0.05]], [7])

        # initial and input sets as sparse polynomial zonotopes
        X0 = convert(SparsePolynomialZonotope, Zonotope(N[1.0, 1.0], N[0.2 0; 0 0.2]))
        U = convert(SparsePolynomialZonotope, Zonotope(zeros(N, 2), N[0.05 0; 0 0.05]))

        # problem definition
        prob = @ivp(x' = A * x + B * u, x(0) ∈ X0, A ∈ AS, B ∈ BS)

        T = N(0.2)
        δ = T / 20

        # recursive true
        alg = HLBS25(; δ=δ,
                     approx_model=CorrectionHullMatrixZonotope(; recursive=true),
                     taylor_order=4,
                     max_order=4,
                     max_order_zono=20,
                     reduction_method=LazySets.GIR05(),
                     recursive=true)

        @test isa(alg, HLBS25)
        @test alg.δ == δ
        @test alg.max_order == 4
        @test alg.max_order_zono == 20
        @test alg.approx_model isa CorrectionHullMatrixZonotope
        @test alg.approx_model.idg === alg.idg

        @test_broken begin
            sol = solve(prob, alg; T=T)
            setrep(sol) <: SparsePolynomialZonotope
        end
        
    end
end
