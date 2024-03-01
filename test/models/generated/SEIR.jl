using ReachabilityAnalysis

@taylorize function seir!(dx, x, p, t)
    S, E, I, R, α, β, γ = x

    βIS = β * (I * S)
    αE = α * E
    γI = γ * I

    dx[1] = -βIS      # dS
    dx[2] = βIS - αE  # dE
    dx[3] = -γI + αE  # dI
    dx[4] = γI        # dR

    # uncertain parameters
    dx[5] = zero(α)
    dx[6] = zero(β)
    dx[7] = zero(γ)
    return dx
end

E₀ = 1e-4
x₀ = [1 - E₀, E₀, 0, 0]
α = 0.2 ± 0.01
β = 1.0 ± 0.0
γ = 0.5 ± 0.01
p = [α, β, γ]
X0 = IntervalBox(vcat(x₀, p));
prob = @ivp(x' = seir!(x), dim:7, x(0) ∈ X0);

sol = solve(prob; T=200.0, alg=TMJets21a(; orderT=7, orderQ=1))
solz = overapproximate(sol, Zonotope);
