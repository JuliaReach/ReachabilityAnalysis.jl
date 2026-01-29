using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector

@testset "LGG09 struct" begin
    # full constructor (requires to manually define the type parameter for the vector type)
    N = Float64
    VN = SingleEntryVector{N}
    TN = BoxDirections{N,VN}
    alg = LGG09{N,Int,VN,TN,Int,Int}(0.01, 0, BoxDirections(1), 0, true, true, true, 0)

    # constructor with default values
    alg = LGG09(; δ=0.01, vars=(1,), dim=1)

    # requires (1) `dirs`, (2) `template`, or (3) `vars` and `dim`
    @test_throws ArgumentError LGG09(; δ=0.01)
    LGG09(; δ=0.01, dirs=[1.0])
    LGG09(; δ=0.01, template=BoxDirections(1))
    @test_throws ArgumentError LGG09(; δ=0.01, vars=(1,))
    @test_throws ArgumentError LGG09(; δ=0.01, dim=1)

    # struct getters
    @test RA.step_size(alg) == 0.01
    @test RA.numtype(alg) == Float64
    @test_throws MethodError setrep(alg)  # TODO fix this
    @test rsetrep(alg) == TemplateReachSet{N,VN,CustomDirections{N,VN},
                                           SubArray{N,1,Matrix{N},Tuple{Base.Slice{Base.OneTo{Int64}},Int64},true}}
end

@testset "LGG09 algorithm: template construction" begin
    # case where the template directions is as single vector
    dirs = fill(0.0, 5)
    dirs[1] = 1.0 # e1
    alg = LGG09(; δ=0.01, template=dirs)
    alg.template == CustomDirections([dirs])

    # case where the template directions is a vector of vectors
    dirs = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(1, 5, -1.0)]
    alg = LGG09(; δ=0.01, template=dirs)
    alg.template == CustomDirections(dirs)

    # alias dirs
    dirs = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(1, 5, -1.0)]
    alg = LGG09(; δ=0.01, dirs=dirs)
    alg.template == CustomDirections(dirs)

    # specification of vars (requires n)
    @test_throws ArgumentError LGG09(δ=0.01, vars=(1))
    alg = LGG09(; δ=0.01, vars=(1), n=5)
    @test collect(alg.template) == dirs

    alg = LGG09(; δ=0.01, vars=(1, 3), n=5)
    dirs_1_3 = [SingleEntryVector(1, 5, 1.0), SingleEntryVector(3, 5, 1.0),
                SingleEntryVector(1, 5, -1.0), SingleEntryVector(3, 5, -1.0)]
    @test collect(alg.template) == dirs_1_3

    # alias using dim
    alg2 = LGG09(; δ=0.01, vars=(1, 3), dim=5)
    @test collect(alg2.template) == dirs_1_3
end

@testset "LGG09 algorithm: 1d" begin
    ivpH, _ = exponential_1d()
    A = ivpH.s.A
    X0 = ivpH.x0
    B = hcat([1.0])
    U = ZeroSet(1)
    X = Universe(1)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)

    δ = 0.01
    box1d = CustomDirections([[1.0], [-1.0]])
    alg = LGG09(; δ=δ, template=box1d)

    for ivp in (ivpH, ivpI)
        # continuous algorithm
        sol = solve(ivp; T=2δ, alg=alg)
        @test isa(sol.alg, LGG09)
        @test setrep(sol) <: HPolyhedron
        @test setrep(sol) == HPolyhedron{Float64,Vector{Float64}}
        @test dim(sol) == 1

        # discrete algorithm
        ivp_norm = ReachabilityAnalysis._normalize(ivp)
        ivp_discr = discretize(ivp_norm, alg.δ, alg.approx_model)
        NSTEPS = 500
        fp_d = ReachabilityAnalysis.post(alg, ivp_discr, NSTEPS)

        # `sparse` option
        alg = LGG09(; δ=δ, template=box1d, sparse=true)
        sol = solve(ivp; T=2δ, alg=alg)
        @test dim(sol) == 1

        # `threaded` option
        alg = LGG09(; δ=δ, template=box1d, threaded=true)
        sol = solve(ivp; T=2δ, alg=alg)
        @test dim(sol) == 1

        # `cache` option deactivated
        alg = LGG09(; δ=δ, template=box1d, cache=false)
        sol = solve(ivp; T=2δ, alg=alg)
        @test dim(sol) == 1
    end

    X = BallInf(zeros(1), 100.0)
    ivpH = @ivp(x' = A * x, x(0) ∈ X0, x ∈ X)
    ivpI = @ivp(x' = A * x + B * u, x(0) ∈ X0, u ∈ U, x ∈ X)

    for ivp in (ivpH, ivpI)
        # invariant
        alg = LGG09(; δ=δ, template=box1d, threaded=false)
        sol = solve(ivp; T=2δ, alg=alg)
        @test dim(sol) == 1
        alg = LGG09(; δ=δ, template=box1d, threaded=true)
        sol = solve(ivp; T=2δ, alg=alg)
        @test dim(sol) == 1
    end

    # `reach_(in)homog_dir_LGG09!` with `Vector` argument (TODO not used internally)
    v = RA.reach_homog_dir_LGG09!([1.0], X0, hcat([1.0]), [1.0], 1, Val{false}())
    @test length(v) == 1
    v = sol = RA.reach_homog_dir_LGG09!([1.0], X0, hcat([1.0]), [1.0], 1, Val{true}())
    @test length(v) == 1
    v = RA.reach_inhomog_dir_LGG09!([1.0], X0, hcat([1.0]), U, [1.0], 1, Val{false}())
    @test length(v) == 1
    v = sol = RA.reach_inhomog_dir_LGG09!([1.0], X0, hcat([1.0]), U, [1.0], 1, Val{true}())
    @test length(v) == 1

    # `reach_homog_eig_LGG09*` (TODO not used internally)
    m = RA.reach_homog_eig_LGG09_posneg([1.0], A, X0, 1)
    @test size(m) == (2, 1)
    m = RA.reach_homog_eig_LGG09([1.0], A, X0, 1)
    @test size(m) == (1, 1)
    t = RA.reach_homog_eig_LGG09_box([1.0], A, X0, 1)
    @test t isa Tuple && size(t[1]) == size(t[2]) == (1, 1)

    # `_upper_bound_eig` (TODO not used internally)
    m = RA._upper_bound_eig(ones(2, 1), hcat([1.0]), 1)
    @test size(m) == (2, 1)
    RA._upper_bound_eig_dir(hcat([1.0]), 1, hcat([1.0]), 1)

    # Krylov
    if isdefined(@__MODULE__, :ExponentialUtilities)
        out = Vector{Float64}(undef, 2)
        RA.reach_homog_krylov_LGG09!(out, X0, A, [1.0], 2)
        RA.reach_inhomog_krylov_LGG09!(out, X0, U, A, [1.0], 2)
    end
end

@testset "LGG09 algorithm: 5d" begin
    # higher-dimensional
    ivp, tspan = linear5D()
    box5d = BoxDirections{Float64,Vector{Float64}}(5)
    sol = solve(ivp; tspan=tspan, alg=LGG09(; δ=0.01, template=box5d))
    @test setrep(sol) == HPolyhedron{Float64,Vector{Float64}}
    @test dim(sol) == 5

    # equivalent algorithm definitions
    alg0 = LGG09(; δ=0.01, template=box5d)
    alg1 = LGG09(; δ=0.01, template=BoxDirections(5))
    alg2 = LGG09(; δ=0.01, template=:box, n=5)
    alg3 = LGG09(; δ=0.01, template=:box, dim=5)
    @test alg1 == alg2 == alg3

    @test !(alg0 == alg1) # alg1-3 use single-entry vectors, but alg0 uses Vector
    @test collect(alg0.template) == Vector.(collect(alg1.template))

    # other directions
    alg4 = LGG09(; δ=0.01, template=:oct, dim=5)
    @test alg4.template isa OctDirections && dim(alg4.template) == 5
end

@testset "LGG09 algorithm: underapproximation" begin
    X0 = BallInf([1.0, 1.0], 0.1)
    function foo()
        A = [0.0 1.0; -1.0 0.0]
        ivp = @ivp(x' = A * x, x(0) ∈ X0)
        tspan = (0.0, 20.0)
        return ivp, tspan
    end
    ivp, tspan = foo()
    δ = 2e-1

    approx_model = SecondOrderddt()
    sol = solve(ivp; tspan=tspan,
                alg=LGG09(; δ=δ, template=PolarDirections(40), approx_model=approx_model))
    approx_model = SecondOrderddt(; oa=false)
    sol_ua = solve(ivp; tspan=tspan,
                   alg=LGG09(; δ=δ, template=PolarDirections(40), approx_model=approx_model))
end
