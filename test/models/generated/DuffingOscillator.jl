using ReachabilityAnalysis

const ω = 1.2

@taylorize function duffing!(du, u, p, t)
    local α = -1.0
    local β = 1.0
    local δ = 0.3
    local γ = 0.37

    x, v = u

    f = γ * cos(ω * t)

    du[1] = v
    du[2] = -α * x - δ * v - β * x^3 + f
    return du
end

X0 = Hyperrectangle([1.0, 0.0], [0.1, 0.1])
prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim:2)

T = 2 * pi / ω;

sol = solve(prob; tspan=(0.0, 20 * T), alg=TMJets21a());
