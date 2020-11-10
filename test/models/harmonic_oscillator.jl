function harmonic_oscillator(; X=Universe(2))
    A = [0.0 1.0; -1.0 0.0]
    prob = @ivp(x' = Ax, x(0) ∈ Singleton([1.0, 0.0]), x ∈ X)
    tspan = (0.0, 20.0)
    return prob, tspan
end
