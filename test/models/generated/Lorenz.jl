using ReachabilityAnalysis

@taylorize function lorenz!(du, u, p, t)
    local σ = 10.0
    local β = 8.0 / 3.0
    local ρ = 28.0

    x, y, z = u
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
    return du
end

X0 = Hyperrectangle(; low=[0.9, 0.0, 0.0], high=[1.1, 0.0, 0.0])
prob = @ivp(x' = lorenz!(x), dim:3, x(0) ∈ X0);

alg = TMJets(; abstol=1e-15, orderT=10, orderQ=2, maxsteps=50_000)
sol = solve(prob; T=10.0, alg=alg)
solz = overapproximate(sol, Zonotope);

@assert ρ([1.0, 0, 0], solz(interval(0, 1.5))) < 20 "the property should be proven"
