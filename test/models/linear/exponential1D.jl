function exponential_1d()
    s = @system(x' = -x)
    x0 = Interval(0.4, 0.5)
    prob = InitialValueProblem(s, x0)
    tspan = (0.0, 1.0)
    return prob, tspan
end
