const A = 1.0
const B = 1.5
const B1 = B + 1

@taylorize function brusselator!(du, u, p, t)
    x, y = u
    x² = x * x
    aux = x² * y
    du[1] = A + aux - B1 * x
    du[2] = B * x - aux
    return du
end

U0(r) = Singleton([1.0, 1.0]) ⊕ BallInf(zeros(2), r)

bruss(r) = @ivp(u' = brusselator!(u), u(0) ∈ U0(r), dim:2)
