@testset "ASB07 algorithm" begin

    # example 1 from
    # "Reachability analysis of linear systems with uncertain parameters and inputs"
    # (Adapted from Reachability.jl)

    # initial set
    X0 = BallInf([1.0, 1.0], 0.1)

    # linear ODE: x' = Ax
    A = IntervalMatrix([-1.0 ± 0.05 -4.0 ± 0.05;
                         4.0 ± 0.05 -1.0 ± 0.05])
    P_lin = @ivp(x' = Ax, x(0) ∈ X0)

    # TODO affine ODE: x' = Ax + Bu
    # B = IntervalMatrix(hcat([1.0 ± 0.01; 1.0 ± 0.0]))
    # u = ConstantInput(Interval(-0.05, 0.05))
    # P_aff = @ivp(x' = Ax + Bu, x(0) ∈ X0)

    f1(P) = solve(P, tspan=(0.0, 1.0), alg=ASB07(δ=0.04));

    f2(P) = solve(P, tspan=(0.0, 5.0), alg=ASB07(δ=0.04, max_order=10),
                  order_discretization = 4);

    f3(P) = solve(P, tspan=(0.0, 1.0), alg=ASB07(δ=0.04),
                  order_discretization = 4, vars =[1, 2],
                  partition = [[1], [2]]);

    f4(P) = solve(P, tspan=(0.0, 5.0), alg=ASB07(δ=0.04, max_order=10),
                  order_discretization = 4, vars =[1, 2],
                  partition = [[1], [2]]);

    for P_t in [P_lin] #[P_lin, P_aff]
        # default options
        for funC in [f1, f2, f3, f4]
            s = funC(P_t)
            @test typeof(s.alg) <: ASB07
            @test setrep(f1(P_lin)) <: Zonotope
        end
    end

end
