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
X0 = Hyperrectangle(vcat(x₀, [0.2, 1, 0.5]), vcat(zeros(4), [0.01, 0, 0.01]))
prob = @ivp(x' = seir!(x), dim:7, x(0) ∈ X0);

sol = solve(prob; T=200.0, alg=TMJets21a(; orderT=7, orderQ=1))
solz = overapproximate(sol, Zonotope);
