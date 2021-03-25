@testset "Reach-set projections" begin
    X = BallInf(zeros(5), 1.0)
    B = BallInf(zeros(2), 1.0)

    R = ReachSet(X, 0 .. 1)

    # concrete projection of a set
    p1 = project(X, [1, 2])
    p2 = project(X, (1, 2))
    p3 = project(X, vars=(1, 2))
    p4 = project(X, vars=[1, 2])
    @test all(x -> isequivalent(B, x), [p1, p2, p3, p4])

    # lazy projection of a set
    p1 = Projection(X, [1, 2])
    p2 = Projection(X, (1, 2))
    p3 = Projection(X, vars=(1, 2))
    p4 = Projection(X, vars=[1, 2])
    @test all(x -> isequivalent(B, overapproximate(x, Hyperrectangle)), [p1, p2, p3, p4])

    # concrete projection of a reach-set
    p1 = project(R, [1, 2])
    p2 = project(R, (1, 2))
    p3 = project(R, vars=(1, 2))
    p4 = project(R, vars=[1, 2])
    @test all(x -> isequivalent(B, set(x)), [p1, p2, p3, p4])

    # concrete projection of a reach-set along a spatial variable and the time variable
    Bt = Hyperrectangle(low=[0.0, -1.0], high=[1.0, 1.0])
    p1 = project(R, [0, 1])
    p2 = project(R, (0, 1))
    p3 = project(R, vars=(0, 1))
    p4 = project(R, vars=[0, 1])
    @test all(x -> isequivalent(Bt, set(x)), [p1, p2, p3, p4])

    # lazy projection of a reach-set
    p1 = Projection(R, [1, 2])
    p2 = Projection(R, (1, 2))
    p3 = Projection(R, vars=(1, 2))
    p4 = Projection(R, vars=[1, 2]);
    @test all(x -> isequivalent(B, overapproximate(set(x), Hyperrectangle)), [p1, p2, p3, p4])

    # lazy projection of a reach-set along a spatial variable and the time variable
    p1 = Projection(R, [0, 1])
    p2 = Projection(R, (0, 1))
    p3 = Projection(R, vars=(0, 1))
    p4 = Projection(R, vars=[0, 1]);
    @test all(x -> isequivalent(Bt, overapproximate(set(x), Hyperrectangle)), [p1, p2, p3, p4])
end

@testset "Conversion of Taylor model reach-sets" begin
    H = Hyperrectangle(ones(2), [0.2, 0.4])
    a = overapproximate(H, TaylorModelReachSet)
    b = convert(TaylorModelReachSet, H)

    @test tspan(a) == tspan(b) == 0 .. 0
    @test isequivalent(set(overapproximate(a, Hyperrectangle)), H)
    @test isequivalent(set(overapproximate(b, Hyperrectangle)), H)

    X = box_approximation(a)
    Y = overapproximate(a, Hyperrectangle)
    @test tspan(X) == tspan(Y)
    @test isequivalent(set(X), set(Y))

    R = ReachSet(H, 0 .. 1)
    c = overapproximate(R, TaylorModelReachSet)
    d = convert(TaylorModelReachSet, R)

    @test tspan(c) == tspan(d) == 0..1
    @test isequivalent(set(overapproximate(c, Hyperrectangle)), set(R))
    @test isequivalent(set(overapproximate(d, Hyperrectangle)), set(R))

    Z = convert(Zonotope, H)
    a = overapproximate(Z, TaylorModelReachSet)
    b = convert(TaylorModelReachSet, Z)

    @test tspan(a) == tspan(b) == 0 .. 0
    @test isequivalent(set(overapproximate(a, Zonotope)), Z)
    @test isequivalent(set(overapproximate(b, Zonotope)), Z)

    R = ReachSet(Z, 0 .. 1)
    c = overapproximate(R, TaylorModelReachSet)
    d = convert(TaylorModelReachSet, R)

    @test tspan(c) == tspan(d) == 0..1
    @test isequivalent(set(overapproximate(c, Zonotope)), set(R))
    @test isequivalent(set(overapproximate(d, Zonotope)), set(R))

end
