using ReachabilityAnalysis, Test
using ReachabilityAnalysis: TimeInterval
import IntervalArithmetic as IA

@testset "TMJets algorithm (TMJets21b)" begin
    prob, tspan = vanderpol()

    # default algorithm for nonlinear systems
    sol = solve(prob; tspan=tspan)
    @test sol.alg isa TMJets

    # pass the algorithm explicitly
    sol = solve(prob; tspan=tspan, alg=TMJets())
    @test sol.alg isa TMJets

    # pass options outside algorithm
    sol = solve(prob; T=0.2, maxsteps=2100)
    @test sol.alg isa TMJets && sol.alg.maxsteps == 2100
    sol = solve(prob; T=0.2, maxsteps=2100) # alias
    @test sol.alg isa TMJets && sol.alg.maxsteps == 2100
    sol = solve(prob; T=0.2, orderT=7, orderQ=1)
    @test sol.alg isa TMJets && sol.alg.orderT == 7 && sol.alg.orderQ == 1
    sol = solve(prob; T=0.2, abstol=1e-11)
    @test sol.alg isa TMJets && sol.alg.abstol == 1e-11
    sol = solve(prob; T=0.2, abstol=1e-11)
    @test sol.alg isa TMJets && sol.alg.abstol == 1e-11

    @test ReachabilityAnalysis.tspan(shift(sol, 1.0)) == ReachabilityAnalysis.tspan(sol) + 1.0

    # split initial conditions
    X0, S = initial_state(prob), system(prob)
    X0s = split(X0, [2, 1]) # split along direction x
    sols = solve(IVP(S, X0s); T=0.1, threading=false) # no threading
    @test flowpipe(sols) isa MixedFlowpipe

    sols = solve(IVP(S, X0s); T=0.1, threading=true) # with threading (default)
    @test flowpipe(sols) isa MixedFlowpipe
end

@ts @testset "TMJets algorithms (TMJets21a/TMJets21b): linear IVPs" begin
    for alg in (TMJets21a, TMJets21b)
        prob, dt = exponential_1d()
        sol = solve(prob; tspan=dt, alg=alg())
        @test sol.alg isa alg

        # getter functions for a taylor model reach-set
        R = sol[1]
        @test TimeInterval(domain(R)) == tspan(R)
        @test IA.diam(remainder(R)[1]) < (alg == TMJets21a ? 1e-13 : 1e-9)
        @test get_order(R) == [8]
        @test polynomial(R) isa Vector{Taylor1{TaylorN{Float64}}}
        @test expansion_point(R) ≈ [IntervalArithmetic.interval(0.0)]

        # test intersection with invariant
        prob, dt = exponential_1d(; invariant=HalfSpace([-1.0], -0.3)) # x >= 0.3
        sol_inv = solve(prob; tspan=dt, alg=alg())
        @test [0.3] ∈ overapproximate(sol_inv[end], Zonotope)
        m = length(sol_inv)
        # check that the following reach-set escapes the invariant
        @test [0.3] ∉ overapproximate(sol[m + 1], Zonotope)
    end

    # TODO test higher order system
    #    prob, tspan = linear5D_homog()
    #    sol = solve(prob, tspan=tspan, alg=TMJets())
    #    @test sol.alg isa TMJets

    # TODO test linear system with input
    #    prob, tspan = linear5D()
    #    sol = solve(prob, tspan=tspan, alg=TMJets())
    #    @test sol.alg isa TMJets
end

@ts @testset "1D Burgers equation (TMJets21b)" begin
    L0 = 1.0 # domain length
    U0 = 1.0 # Re = 20.
    x = range(-0.5 * L0, 0.5 * L0; length=4)
    # Initial velocity
    X0 = Singleton(-U0 * sin.(2 * π / L0 * x))
    # IVP definition
    prob = @ivp(x' = burgers!(x), dim = 4, x(0) ∈ X0)
    sol = solve(prob; tspan=(0.0, 1.0), alg=TMJets21b())
    @test dim(sol) == 4
end

#=
alg = TMJets(abstol=1e-10, orderT=10, orderQ=2)

# reach mode
sol = solve(P, T=7.0, alg=alg)
@test set(sol[1]) isa Hyperrectangle # check default set representation
=#
