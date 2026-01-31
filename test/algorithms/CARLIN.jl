using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
import IntervalArithmetic as IA
import Symbolics

@testset "CARLIN struct" begin
    # full constructor
    alg = CARLIN(1, true, 0.01, true, 0)

    # constructor with default values
    alg = CARLIN(; δ=0.01)

    # struct getters
    @test_broken RA.numtype(alg) == Float64  # TODO this should be supported
    @test_broken setrep(alg) <: LazySet   # TODO this should be supported
    @test_broken rsetrep(alg) <: ReachSet  # TODO this should be supported
end

@testset "CARLIN lift_vector" begin
    x = IA.interval(0, 1)
    H = RA.lift_vector(x, 2)
    H2 = RA.lift_vector(Interval(x), 2)
    @test H == H2 == Hyperrectangle([0.5, 0.5], [0.5, 0.5])
end

@testset "CARLIN kron_pow" begin
    # Interval
    X0 = Interval(0.47, 0.53)
    for N in (1, 2)
        I1 = RA.kron_pow(X0, N)
        I2 = RA.kron_pow(X0.dat, N)  # IA.Interval
        @test I1 == Interval(I2)

        if N == 1
            @test RA.kron_pow(X0, N, "invalid") == X0
        else
            @test_throws ArgumentError RA.kron_pow(X0, N, "invalid")
        end

        # TODO "explicit" algorithm cannot be passed to `Interval` method,
        # so it chooses the `Hyperrectangle` method
        I2 = RA.kron_pow(X0, N, "explicit")
        if N == 2
            @test_broken I2 isa Interval
        end
        @test isequivalent(I2, I1)
    end

    # Hyperrectangle
    X0 = Hyperrectangle([1.0, 1, 1], [0.1, 0.1, 0.1])
    # explicit `kron_pow` (TODO refactor code so that this can be chosen by the user)
    H = RA.kron_pow(X0, 2, "explicit")
    @test H ≈ Hyperrectangle(fill(1.01, 9), fill(0.2, 9))
end

@testset "CARLIN kron_pow_stack" begin
    # Interval
    X0 = Interval(0.47, 0.53)
    N = 2
    H1 = RA.kron_pow_stack(X0, N)
    H2 = RA.kron_pow_stack(X0.dat, N)  # IA.Interval
    @test H1 == H2

    # Hyperrectangle
    X0 = Hyperrectangle([1.0, 1, 1], [0.1, 0.1, 0.1])
    H = RA.kron_pow_stack(X0, N)
    @test_broken H isa Hyperrectangle  # TODO `kron_pow_stack` is not consistent in return type
end

@testset "CARLIN algorithm: logistic model" begin
    r = -0.5
    K = 0.8
    F1 = hcat(r)
    F2 = hcat(-r / K)
    X0 = Interval(0.47, 0.53)
    prob = @ivp(CanonicalQuadraticForm(F1, F2), x(0) ∈ X0)

    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=true))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.2

    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=false))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.01

    # Use the compressed Kronecker form
    sol = solve(prob; T=10.0, alg=CARLIN(; N=3, bloat=false, compress=true))
    @test dim(sol) == 1
    @test ρ([1.0], sol(10.0)) < 0.01
end

@testset "CARLIN algorithm: SEIR model" begin
    F1 = zeros(3, 3)
    F1[1, 1] = -0.19
    F1[2, 2] = -0.192308
    F1[3, 2] = 0.192308
    F1[3, 3] = -0.434783
    F2 = zeros(3, 9)
    F2[1, 3] = -1.3e-8
    F2[2, 3] = 1.3e-8
    X0 = Hyperrectangle([1.0, 1, 1], [0.1, 0.1, 0.1])
    prob = @ivp(CanonicalQuadraticForm(F1, F2), x(0) ∈ X0)

    sol = solve(prob; T=1.0, alg=CARLIN(; N=2, bloat=true))
    @test dim(sol) == 3

    # reset option: uniform resets
    sol = solve(prob; T=1.0, alg=CARLIN(; resets=2))
    @test dim(sol) == 3

    # reset option: specific resets
    sol = solve(prob; T=1.0, alg=CARLIN(; resets=[0.1]))
    @test dim(sol) == 3
end

@testset "CARLIN algorithm: black-box model" begin
    function dyn!(dx, x, p, t)
        dx[1] = -x[1]
        return dx
    end
    X0 = Interval(0.4, 0.5)
    prob = @ivp(x' = dyn!(x), dim:1, x(0) ∈ X0)
    @test_broken solve(prob; T=10.0, alg=CARLIN(; N=3)) isa Flowpipe  # TODO finish implementation
end
