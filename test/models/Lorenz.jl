using ReachabilityAnalysis, Plots

@taylorize function lorenz!(dx, x, params, t)
    local σ = 10.0
    local β = 8.0 / 3.0
    local ρ = 28.0
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
    return dx
end

X0 = Hyperrectangle(low=[0.9, 0.0, 0.0], high=[1.1, 0.0, 0.0])
prob = @ivp(x' = lorenz!(x), dim=3, x(0) ∈ X0);

alg = TMJets(abstol=1e-15, orderT=10, orderQ=2, maxsteps=50_000);

sol = solve(prob, T=10.0, alg=alg);

solz = overapproximate(sol, Zonotope);

plot(solz, vars=(0, 1), xlab="t", ylab="x")

plot(solz(0.0 .. 1.5), vars=(0, 1), xlab="t", ylab="x", lw=0.0)
plot!(x -> 20.0, c=:red, xlims=(0.0, 1.5), lab="")

ρ([1.0, 0.0, 0.0], solz(0.0 .. 1.5))

plot(solz, vars=(0, 2), xlab="t", ylab="y")

-ρ([0.0, -1.0, 0.0], solz)

ρ([0.0, 1.0, 0.0], solz)

plot(solz, vars=(0, 3), xlab="t", ylab="z")

plot(solz, vars=(1, 3))

