using ReachabilityAnalysis, Test
import ReachabilityAnalysis as RA
import Flowstar, Symbolics
using Suppressor: @suppress  # suppress long output of Flowstar

@testset "FLOWSTAR struct" begin
    # full constructor
    alg = FLOWSTAR(nothing, nothing, 0.01, nothing, 0.0, 0, true, nothing)

    # constructor with default values
    alg = FLOWSTAR(; δ=0.01)

    # struct getters
    @test RA.numtype(alg) == Float64
    @test_broken setrep(alg) == Vector{TaylorModel1{TaylorN{Float64},Float64}}  # TODO this should work
    @test rsetrep(alg) == TaylorModelReachSet{Float64,RA.IA.Interval{Float64}}
end

@testset "FLOWSTAR algorithm" begin
    # Passing a model file.
    model = joinpath(@__DIR__, "..", "models", "LotkaVolterra.model")
    ivp = @ivp(BlackBoxContinuousSystem(model, 2), x(0) ∈ Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2]))
    sol = @suppress solve(ivp; tspan=(0, 1), alg=FLOWSTAR(; δ=0.02))
    RT = TaylorModelReachSet{Float64,IntervalArithmetic.Interval{Float64}}
    @test sol isa RA.ReachSolution{Flowpipe{Float64,RT,Vector{RT}},
                                   FLOWSTAR{Float64,Flowstar.FixedTMOrder,
                                            Flowstar.QRPreconditioner,
                                            Flowstar.NonPolyODEScheme}}

    # Generating the model file from the given system.
    function f!(dx, x, p, t)
        dx[1] = 1.5 * x[1] - x[1] * x[2]
        return dx[2] = -3 * x[2] + x[1] * x[2]
    end
    ivp = @ivp(x' = f!(x), dim = 2, x(0) ∈ Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2]))
    sol = @suppress solve(ivp; tspan=(0, 1), alg=FLOWSTAR(; δ=0.02))
    RT = TaylorModelReachSet{Float64,IntervalArithmetic.Interval{Float64}}
    @test sol isa RA.ReachSolution{Flowpipe{Float64,RT,Vector{RT}},
                                   FLOWSTAR{Float64,Flowstar.FixedTMOrder,
                                            Flowstar.QRPreconditioner,
                                            Flowstar.NonPolyODEScheme}}
end
