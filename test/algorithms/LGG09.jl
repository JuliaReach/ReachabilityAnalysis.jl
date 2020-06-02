using LazySets.Arrays: SingleEntryVector

@testset "LGG09 algorithm" begin
    # one-dimensional
    prob, tspan = exponential_1d()
    box1d = CustomDirections([[1.0], [-1.0]])
    sol = solve(prob, tspan=tspan, LGG09(δ=0.01, template=box1d))
    @test isa(sol.alg, LGG09)
    @test setrep(sol) <: HPolyhedron
    @test setrep(sol) == HPolyhedron{Float64, Vector{Float64}}
    @test dim(sol) == 1

    # higher-dimensional
    prob, tspan = linear5D()
    box5d = BoxDirections{Float64, Vector{Float64}}(5)
    sol = solve(prob, tspan=tspan, LGG09(δ=0.01, template=box5d))
    @test setrep(sol) == HPolyhedron{Float64, Vector{Float64}}
    @test dim(sol) == 5

    # case where the template directions is as single vector
    dirs = fill(0.0, 5); dirs[1] = 1.0 # e1
    alg = LGG09(δ=0.01, template=dirs)
    alg.template == CustomDirections([dirs])

    # case where the template directions is a vector of vectors
    dirs = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(1, 5, -1.0)]
    alg = LGG09(δ=0.01, template=dirs)
    alg.template == CustomDirections(dirs)
end
