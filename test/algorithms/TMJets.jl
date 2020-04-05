@testset "TMJets algorithm" begin
    prob, tspan = vanderpol()

    # default algorithm for nonlinear systems
    sol = solve(prob, tspan=tspan)
    @test sol.alg isa TMJets

    # pass the algorithm explicitly
    sol = solve(prob, tspan=tspan, TMJets())
    @test sol.alg isa TMJets

    # TODO: try different options
end

@testset "TMJets algorithm: linear IVPs" begin
    prob, tspan = exponential_1d()
    sol = solve(prob, tspan=tspan, TMJets())
    @test sol.alg isa TMJets

    # TODO test higher order system
#    prob, tspan = linear5D_homog()
#    sol = solve(prob, tspan=tspan, TMJets())
#    @test sol.alg isa TMJets

    # TODO test linear system with input
#    prob, tspan = linear5D()
#    sol = solve(prob, tspan=tspan, TMJets())
#    @test sol.alg isa TMJets
end

#=
alg = TMJets(abs_tol=1e-10, orderT=10, orderQ=2)

# reach mode
sol = solve(P, T=7.0, alg)
@test set(sol[1]) isa Hyperrectangle # check default set representation
=#
