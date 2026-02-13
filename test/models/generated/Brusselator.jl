using ReachabilityAnalysis
import ReachabilityAnalysis.ReachabilityBase.Comparison as CMP

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

U₀ = Hyperrectangle(low=[0.8, 0], high=[1, 0.2])
prob = @ivp(u' = brusselator!(u), u(0) ∈ U₀, dim:2)

T = 18.0;

alg = TMJets21a(; orderT=6, orderQ=2)
sol = solve(prob; T=T, alg=alg);

U0(r) = BallInf([1.0, 1.0], r);

bruss(r) = @ivp(u' = brusselator!(u), u(0) ∈ U0(r), dim:2);

sol_01 = solve(bruss(0.01); T=30.0, alg=alg);

sol_1 = solve(bruss(0.1); T=30.0, alg=alg);
