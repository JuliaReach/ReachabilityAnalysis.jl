using ReachabilityAnalysis: _jacobian_function,
                            _jacobian_expression,
                            _hessian_function,
                            _hessian_expression

@testset "Symbolic-Numeric utils: Jacobian" begin
    # Brusselator model as a vector of symbolic expressions
    var = @variables x y
    f = [1 + x^2 * y - 2.5x, 1.5x - x^2 * y]
    Jexpr0 = [-2.5+2x * y x^2; 1.5-2x * y -x^2]

    Jexpr = _jacobian_expression(f, var)
    @test isequal(Jexpr, Jexpr0)
    J = _jacobian_function(f, var)
    v = [1.0, 2.0]
    Jv = [1.5 1; -2.5 -1]
    @test J(v) == Jv

    # same using static arrays
    f_st = SVector{2,Num}(f)
    Jexpr_st = _jacobian_function(f_st, var)
    J_st = _jacobian_function(f_st, var)
    v_st = SA[1.0, 2.0]
    Jv_st = SMatrix{2,2,Float64}(Jv)
    aux = J_st(v_st)
    @test isa(aux, StaticArray) && aux == Jv_st

    # same using MathematicalSystem types
    sys = @system(u' = f(u), dim = 2)
    J = _jacobian_function(sys, var)
    @test J(v) == Jv
    Jexpr = _jacobian_expression(sys, var)
    @test isequal(Jexpr0, Jexpr)

    U0 = VPolygon([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    ivp = @ivp(u' = f(u), dim = 2, u(0) ∈ U0)
    J = _jacobian_function(ivp, var)
    @test J(v) == Jv
    Jexpr = _jacobian_expression(ivp, var)
    @test isequal(Jexpr0, Jexpr)

    # test in-place version
    prob_iip = bruss(0.1)
    J = _jacobian_function(prob_iip, var)
    @test J(v) == Jv
    Jexpr = _jacobian_expression(prob_iip, var)
    @test isequal(Jexpr0, Jexpr)

    # linear systems
    A = [1 3; 2 -2.0]
    sys = @system(x' = Ax)
    J = _jacobian_function(sys, var)
    aux = [1.0, 1.0] # can be any vector
    @test J(aux) == A
    Jexpr = _jacobian_expression(sys, var)
    @test isequal(Jexpr, A)
end

@testset "Symbolic-Numeric utils: Hessian" begin
end
