function harmonic_oscillator(; X=Universe(2))
    A = [0.0 1.0; -1.0 0.0]
    prob = @ivp(x' = Ax, x(0) ∈ Singleton([1.0, 0.0]), x ∈ X)
    tspan = (0.0, 20.0)
    return prob, tspan
end

function harmonic_oscillator_homog_RFEM()
    A = [0 1; -(4π)^2 0]
    X0 = Hyperrectangle([1.0, 0.0], [0.1, 0.1])
    prob = @ivp(x' = A*x, x(0) ∈ X0)
    tspan = (0.0, 2.0)
    return prob, tspan
end

function harmonic_oscillator_forced_RFEM()
    A = [0 1; -(4π)^2 0]
    X0 = Hyperrectangle([1.0, 0.0], [0.1, 0.1])
    B = Matrix(1.0I, 2, 2)
    X = Universe(2)
    U = Interval(0.0) × Interval(0.8, 1.2)
    prob = @ivp(x' = A*x + B*u, x(0) ∈ X0, x ∈ X, u ∈ U)
    tspan = (0.0, 2.0)
    return prob, tspan
end
