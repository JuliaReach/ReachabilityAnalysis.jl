using ReachabilityAnalysis

@taylorize function laubloomis!(dx, x, p, t)
    dx[1] = 1.4 * x[3] - 0.9 * x[1]
    dx[2] = 2.5 * x[5] - 1.5 * x[2]
    dx[3] = 0.6 * x[7] - 0.8 * (x[2] * x[3])
    dx[4] = 2 - 1.3 * (x[3] * x[4])
    dx[5] = 0.7 * x[1] - (x[4] * x[5])
    dx[6] = 0.3 * x[1] - 3.1 * x[6]
    dx[7] = 1.8 * x[6] - 1.6 * (x[2] * x[7])
    return dx
end

function laubloomis(; W=0.01)
    X0 = Hyperrectangle([1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45], fill(W, 7))
    prob = @ivp(x' = laubloomis!(x), dim:7, x(0) ∈ X0)
    return prob
end;

const e4 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

prob = laubloomis(; W=0.01)
alg = TMJets(; abstol=1e-11, orderT=7, orderQ=1, adaptive=true)

sol_1 = solve(prob; T=20.0, alg=alg)
sol_1z = overapproximate(sol_1, Zonotope);

@assert ρ(e4, sol_1z) < 4.5 "the property should be proven"

prob = laubloomis(; W=0.05)
alg = TMJets(; abstol=1e-12, orderT=7, orderQ=1, adaptive=false)

sol_2 = solve(prob; T=20.0, alg=alg)
sol_2z = overapproximate(sol_2, Zonotope);

@assert ρ(e4, sol_2z) < 4.5 "the property should be proven"

prob = laubloomis(; W=0.1)
alg = TMJets(; abstol=1e-12, orderT=7, orderQ=1, adaptive=false)

sol_3 = solve(prob; T=20.0, alg=alg)
sol_3z = overapproximate(sol_3, Zonotope);

@assert ρ(e4, sol_3z) < 5.0 "the property should be proven"
