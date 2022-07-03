using LazySets.Arrays: SingleEntryVector

@testset "LGG09 algorithm: template construction" begin
    # case where the template directions is as single vector
    dirs = fill(0.0, 5); dirs[1] = 1.0 # e1
    alg = LGG09(δ=0.01, template=dirs)
    alg.template == CustomDirections([dirs])

    # case where the template directions is a vector of vectors
    dirs = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(1, 5, -1.0)]
    alg = LGG09(δ=0.01, template=dirs)
    alg.template == CustomDirections(dirs)

    # alias dirs
    dirs = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(1, 5, -1.0)]
    alg = LGG09(δ=0.01, dirs=dirs)
    alg.template == CustomDirections(dirs)

    # specification of vars (requires n)
    @test_throws ArgumentError LGG09(δ=0.01, vars=(1))
    alg = LGG09(δ=0.01, vars=(1), n=5)
    @test collect(alg.template) == dirs

    alg = LGG09(δ=0.01, vars=(1, 3), n=5)
    dirs_1_3 = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(3, 5, 1.0),
                SingleEntryVector(1, 5, -1.0), SingleEntryVector(3, 5, -1.0)]
    @test collect(alg.template) == dirs_1_3

    # alias using dim
    alg2 = LGG09(δ=0.01, vars=(1, 3), dim=5)
    @test collect(alg2.template) == dirs_1_3
end

@testset "LGG09 algorithm: 1d" begin
    # one-dimensional
    prob, dt = exponential_1d()
    box1d = CustomDirections([[1.0], [-1.0]])
    sol = solve(prob, tspan=dt, LGG09(δ=0.01, template=box1d))
    @test isa(sol.alg, LGG09)
    @test setrep(sol) <: HPolyhedron
    @test setrep(sol) == HPolyhedron{Float64, Vector{Float64}}
    @test dim(sol) == 1
end

@testset "LGG09 algorithm: 5d" begin
    # higher-dimensional
    prob, dt = linear5D()
    box5d = BoxDirections{Float64, Vector{Float64}}(5)
    sol = solve(prob, tspan=dt, LGG09(δ=0.01, template=box5d))
    @test setrep(sol) == HPolyhedron{Float64, Vector{Float64}}
    @test dim(sol) == 5

    # equivalent algorithm definitions
    alg0 = LGG09(δ=0.01, template=box5d)
    alg1 = LGG09(δ=0.01, template=BoxDirections(5))
    alg2 = LGG09(δ=0.01, template=:box, n=5)
    alg3 = LGG09(δ=0.01, template=:box, dim=5)
    @test alg1 == alg2 == alg3

    @test !(alg0 == alg1) # alg1-3 use single-entry vectors, but alg0 uses Vector
    @test collect(alg0.template) == Vector.(collect(alg1.template))
    
    # other directions
    alg4 = LGG09(δ=0.01, template=:oct, dim=5)
    @test alg4.template isa OctDirections && dim(alg4.template) == 5
end

@testset "LGG09 algorithm: underapproximation" begin
    X0 = BallInf([1.0, 1.0], 0.1)
    function foo()
        A = [0.0 1.0; -1.0 0.0]
        prob = @ivp(x' = Ax, x(0) ∈ X0)
        tspan = (0.0, 20.0)
        return prob, tspan
    end
    prob, dt = foo()
    δ = 2e-1

    approx_model = SecondOrderddt()
    sol = solve(prob, tspan=(0.0, 20.0),
                alg=LGG09(δ=δ, template=PolarDirections(40), approx_model=approx_model))
    approx_model = ReachabilityAnalysis.SecondOrderddt(oa=false)
    sol_ua = solve(prob, tspan=(0.0, 20.0),
                   alg=LGG09(δ=δ, template=PolarDirections(40), approx_model=approx_model))
end
