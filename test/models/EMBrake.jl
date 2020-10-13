set_rtol(Float64, 1e-12)
set_ztol(Float64, 1e-13)

# Fixed parameters

function embrake_no_pv(; Tsample=1.E-4, ζ=1e-6, x0=0.05)
    # model's constants
    L = 1.e-3
    KP = 10000.
    KI = 1000.
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167

    # state variables: [I, x, xe, xc]
    A = Matrix([-(R+K^2/drot)/L 0      KP/L  KI/L;
                          K/i/drot      0      0      0;
                          0             0      0      0;
                          0             0      0      0])

    EMbrake = @system(x' = Ax)

    # initial conditions
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = I₀ × x₀ × xe₀ × xc₀

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(EMbrake, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

# Parameter variation

function embrake_pv_1(; Tsample=1.E-4, ζ=1e-6, Δ=3.0, x0=0.05)
    # model's constants
    L = 1.e-3
    KP = 10000.
    KI = 1000.
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167
    p = 504. + (-Δ .. Δ)

    # state variables: [I, x, xe, xc]
    A = IntervalMatrix([-p            0      KP/L   KI/L;
                        K/i/drot      0      0      0;
                        0             0      0      0;
                        0             0      0      0])

    EMbrake = @system(x' = Ax)

    # initial conditions
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = I₀ × x₀ × xe₀ × xc₀

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(EMbrake, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

# Extended parameter variation

function embrake_pv_2(; Tsample=1.E-4, ζ=1e-6, x0=0.05, χ=5.0)
    # model's constants
    Δ = -χ/100 .. χ/100
    L = 1.e-3 * (1 + Δ)
    KP = 10000. * (1 + Δ)
    KI = 1000. * (1 + Δ)
    R = 0.5 * (1 + Δ)
    K = 0.02 * (1 + Δ)
    drot = 0.1 * (1 + Δ)
    i = 113.1167 * (1 + Δ)

    # state variables: [I, x, xe, xc]
    A = IntervalMatrix([-(R+K^2/drot)/L 0      KP/L  KI/L;
                          K/i/drot      0      0      0;
                          0             0      0      0;
                          0             0      0      0])

    EMbrake = @system(x' = Ax)

    # initial conditions
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = I₀ × x₀ × xe₀ × xc₀

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(EMbrake, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

