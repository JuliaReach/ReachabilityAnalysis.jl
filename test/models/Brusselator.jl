using ReachabilityAnalysis

const A = 1.0
const B = 1.5
const B1 = B + 1

@taylorize function brusselator!(du, u, p, t)
    x, y = u
    x² = x * x
    aux = x² * y
    du[1] = A + aux - B1*x
    du[2] = B*x - aux
    return du
end

U0(r) = Singleton([1.0, 1.0]) ⊕ BallInf(zeros(2), r)

bruss(r) = @ivp(u' = brusselator!(u), u(0) ∈ U0(r), dim: 2)

sol_01 = solve(bruss(0.01), T=30.0, alg=TMJets(orderT=6, orderQ=2))

LazySets.set_ztol(Float64, 1e-15)

sol_1 = solve(bruss(0.1), T=30.0, alg=TMJets(orderT=6, orderQ=2))

vol_01 = overapproximate(sol_01(9.0), Hyperrectangle) |> set |> volume

vol_1 = overapproximate(sol_1(9.0), Hyperrectangle) |> set |> volume

