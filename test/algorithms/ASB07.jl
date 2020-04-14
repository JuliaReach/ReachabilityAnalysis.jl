@testset "ASB07 algorithm" begin
    # example 1 from
    # "Reachability analysis of linear systems with uncertain parameters and inputs"

    # initial set
    X0 = BallInf([1.0, 1.0], 0.1)

    # linear ODE: x' = Ax
    A = IntervalMatrix([-1.0 ± 0.05 -4.0 ± 0.05;
                         4.0 ± 0.05 -1.0 ± 0.05])
    P_lin = @ivp(x' = Ax, x(0) ∈ X0)

    sol1 = solve(P_lin, tspan=(0.0, 1.0), alg=ASB07(δ=0.04));
    @test dim(sol1)==2
    @test isa(sol1.alg, ASB07)
    @test sol1.alg.δ==0.04
    @test sol1.alg.max_order==10
    @test typeof(sol1.alg.approx_model)==CorrectionHull
    @test sol1.alg.approx_model.order==10
    @test setrep(sol1) <: Zonotope
    @test setrep(sol1)==Zonotope{Float64,Array{Float64,1},Array{Float64,2}}

    sol2 = solve(P_lin, tspan=(0.0, 5.0), alg=ASB07(δ=0.01, max_order=7,
                 approx_model=CorrectionHull(order=9, exp=:base)));
    @test dim(sol2)==2
    @test isa(sol2.alg, ASB07)
    @test sol2.alg.δ==0.01
    @test sol2.alg.max_order==7
    @test typeof(sol2.alg.approx_model)==CorrectionHull
    @test sol2.alg.approx_model.order==9
    @test setrep(sol2) <: Zonotope
    @test setrep(sol2)==Zonotope{Float64,Array{Float64,1},Array{Float64,2}}

end

# TODO affine ODE: x' = Ax + Bu
# B = IntervalMatrix(hcat([1.0 ± 0.01; 1.0 ± 0.0]))
# u = ConstantInput(Interval(-0.05, 0.05))
# P_aff = @ivp(x' = Ax + Bu, x(0) ∈ X0)
