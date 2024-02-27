@taylorize function brusselator!(du, u, p, t)
    local A = 1.0
    local B = 1.5
    local B1 = B + 1

    x, y = u

    x²y = x^2 * y
    du[1] = A + x²y - B1 * x
    du[2] = B * x - x²y
    return du
end

U0(r) = BallInf([1.0, 1.0], r);

bruss(r) = @ivp(u' = brusselator!(u), u(0) ∈ U0(r), dim:2);
